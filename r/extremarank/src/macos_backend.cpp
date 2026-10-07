#include "count_reader.h"
#include <algorithm>
#include <cstring>
#if defined(__APPLE__)
#import <Foundation/Foundation.h>
#import <Metal/Metal.h>
#endif
using namespace Rcpp;

#if defined(__APPLE__)
static id<MTLDevice> er_device() {
    static id<MTLDevice> device=MTLCreateSystemDefaultDevice(); return device;
}
static id<MTLComputePipelineState> er_pipeline(id<MTLDevice>d) {
    static id<MTLComputePipelineState> pipeline=nil;
    if(pipeline) return pipeline;
    // Every count is split into two integer limbs. Carry propagation uses the
    // returned old low limb, so concurrent additions remain exact.
    NSString *source=@R"METAL(
#include <metal_stdlib>
using namespace metal;
struct Entry { uint gene, unit, lo, hi; };
struct Sum { atomic_uint lo; atomic_uint hi; };
kernel void grouped_sum(device const Entry *entries [[buffer(0)]],
    device Sum *sums [[buffer(1)]], device atomic_uint *overflow [[buffer(2)]],
    constant uint &n [[buffer(3)]], constant uint &genes [[buffer(4)]],
    uint i [[thread_position_in_grid]]) {
    if(i>=n) return;
    Entry e=entries[i]; uint target=e.unit*genes+e.gene;
    uint old=atomic_fetch_add_explicit(&sums[target].lo,e.lo,memory_order_relaxed);
    uint high=e.hi+uint(old>0xffffffffu-e.lo);
    if(high) {
        uint old_high=atomic_fetch_add_explicit(&sums[target].hi,high,memory_order_relaxed);
        if(old_high>0xffffffffu-high) atomic_store_explicit(overflow,1u,memory_order_relaxed);
    }
}
)METAL";
    NSError *error=nil;
    id<MTLLibrary> lib=[d newLibraryWithSource:source options:nil error:&error];
    if(!lib) stop("Metal shader compilation failed: %s",[[error localizedDescription]UTF8String]);
    id<MTLFunction> f=[lib newFunctionWithName:@"grouped_sum"];
    pipeline=[d newComputePipelineStateWithFunction:f error:&error];
    if(!pipeline) stop("Metal compute pipeline failed: %s",[[error localizedDescription]UTF8String]);
    return pipeline;
}
#endif

// [[Rcpp::export]]
List cpp_backend_info() {
#if defined(__APPLE__)
    @autoreleasepool {
        id<MTLDevice>d=er_device();
        if(!d) return List::create(_["metal_available"]=false,_["reason"]="Metal device unavailable (including restricted execution contexts)",_["native_available"]=true,_["cuda_available"]=false);
        return List::create(_["metal_available"]=true,_["device"]=std::string([[d name]UTF8String]),
            _["unified_memory"]=(bool)[d hasUnifiedMemory],_["max_buffer_bytes"]=(double)[d maxBufferLength],
            _["recommended_working_set_bytes"]=(double)[d recommendedMaxWorkingSetSize],
            _["native_available"]=true,_["cuda_available"]=false);
    }
#else
    return List::create(_["metal_available"]=false,_["reason"]="Metal backend requires macOS",_["native_available"]=true,_["cuda_available"]=false);
#endif
}

// [[Rcpp::export]]
List cpp_pseudobulk_metal(SEXP x,IntegerVector unit,int n_units,double max_gpu_bytes) {
#if defined(__APPLE__)
    @autoreleasepool {
        id<MTLDevice>d=er_device(); if(!d) stop("Metal device unavailable; use backend='native'");
        ercounts::Reader reader(x); ercounts::check_units(reader,unit,n_units);
        size_t outputs=static_cast<size_t>(reader.nr)*n_units;
        if(outputs>std::numeric_limits<uint32_t>::max()) stop("Metal output indexing limit exceeded; use backend='native'");
        size_t output_bytes=outputs*8, fixed=output_bytes+4;
        if(!std::isfinite(max_gpu_bytes) || max_gpu_bytes<static_cast<double>(fixed+sizeof(ercounts::Entry)*2))
            stop("max_gpu_bytes must cover integer outputs and at least two input entries");
        if(output_bytes>[d maxBufferLength]) stop("Metal output exceeds device buffer limit");
        double remaining=max_gpu_bytes-fixed;
        size_t entries=std::max<R_xlen_t>(reader.length,1);
        int slots=entries<=remaining/sizeof(ercounts::Entry)?1:2;
        size_t capacity=std::min<size_t>(entries,std::min<double>(remaining/(slots*sizeof(ercounts::Entry)),
            std::min<double>([d maxBufferLength]/sizeof(ercounts::Entry),std::numeric_limits<uint32_t>::max())));
        if(!capacity) stop("Metal staging budget too small");
        id<MTLBuffer> out=[d newBufferWithLength:output_bytes options:MTLResourceStorageModeShared];
        id<MTLBuffer> flag=[d newBufferWithLength:4 options:MTLResourceStorageModeShared];
        id<MTLBuffer> inputs[2]={nil,nil};
        id<MTLCommandBuffer> pending[2]={nil,nil};
        if(!out || !flag) stop("Metal output allocation failed");
        std::memset([out contents],0,output_bytes); std::memset([flag contents],0,4);
        for(int s=0;s<slots;++s) { inputs[s]=[d newBufferWithLength:capacity*sizeof(ercounts::Entry) options:MTLResourceStorageModeShared]; if(!inputs[s]) stop("Metal staging allocation failed"); }
        id<MTLComputePipelineState> pipeline=er_pipeline(d);
        id<MTLCommandQueue> queue=[d newCommandQueue]; if(!queue) stop("Metal queue creation failed");
        double gpu_seconds=0; int commands=0; R_xlen_t scanned=0; bool finished=false;
        auto wait=[&](int s){
            if(pending[s]) {
                [pending[s]waitUntilCompleted];
                if([pending[s]status]==MTLCommandBufferStatusError) stop("Metal execution failed: %s",[[[pending[s]error]localizedDescription]UTF8String]);
                double begin=[pending[s]GPUStartTime],end=[pending[s]GPUEndTime]; if(end>=begin && begin>0) gpu_seconds+=end-begin;
                pending[s]=nil;
            }
        };
        while(!finished) {
            int s=commands%slots; wait(s); checkUserInterrupt();
            auto *buffer=static_cast<ercounts::Entry *>([inputs[s]contents]); size_t n=0;
            int g,c; double value;
            while(n<capacity) {
                if(!reader.next(g,c,value)) { finished=true; break; }
                if((++scanned & 1048575)==0) checkUserInterrupt();
                int u=unit[c]; if(!u || !value) continue;
                uint64_t v=static_cast<uint64_t>(value);
                buffer[n++]={static_cast<uint32_t>(g),static_cast<uint32_t>(u-1),static_cast<uint32_t>(v),static_cast<uint32_t>(v>>32)};
            }
            if(!n) break;
            id<MTLCommandBuffer> command=[queue commandBuffer];
            id<MTLComputeCommandEncoder> encoder=[command computeCommandEncoder];
            if(!command || !encoder) stop("Metal command creation failed");
            [encoder setComputePipelineState:pipeline];
            [encoder setBuffer:inputs[s] offset:0 atIndex:0]; [encoder setBuffer:out offset:0 atIndex:1]; [encoder setBuffer:flag offset:0 atIndex:2];
            uint32_t count=n,genes=reader.nr;
            [encoder setBytes:&count length:4 atIndex:3]; [encoder setBytes:&genes length:4 atIndex:4];
            NSUInteger threads=std::min<NSUInteger>(256,[pipeline maxTotalThreadsPerThreadgroup]);
            [encoder dispatchThreadgroups:MTLSizeMake((n+threads-1)/threads,1,1) threadsPerThreadgroup:MTLSizeMake(threads,1,1)];
            [encoder endEncoding]; [command commit]; pending[s]=command; ++commands;
        }
        for(int s=0;s<slots;++s) wait(s);
        if(*static_cast<uint32_t *>([flag contents])) stop("pseudobulk integer accumulator overflow; sums exceed exact range");
        auto *values=static_cast<uint32_t *>([out contents]); NumericMatrix result(reader.nr,n_units);
        for(size_t i=0;i<outputs;++i) {
            uint64_t value=(static_cast<uint64_t>(values[2*i+1])<<32)|values[2*i];
            if(value>static_cast<uint64_t>(ercounts::max_exact)) stop("pseudobulk sum exceeds exact integer range 2^53 - 1");
            result[i]=static_cast<double>(value);
        }
        return List::create(_["counts"]=result,_["backend"]="metal",_["device"]=std::string([[d name]UTF8String]),
            _["gpu_buffer_bytes"]=static_cast<double>(fixed+slots*capacity*sizeof(ercounts::Entry)),
            _["gpu_command_seconds"]=gpu_seconds,_["gpu_commands"]=commands,
            _["stored_values_scanned"]=static_cast<double>(scanned),_["unified_memory"]=(bool)[d hasUnifiedMemory]);
    }
#else
    stop("Metal backend requires macOS; use backend='native'"); return List();
#endif
}
