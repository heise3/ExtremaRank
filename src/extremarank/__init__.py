"""ExtremaRank: exact studentized subset extrema and shared-donor rank audit."""
from .core import Score, ExtremaResult, PreparedValues, compare_scores, compare_abs_scores, compare_t, extrema, score
from .audit import AuditResult, audit_topk

__version__ = "0.1.0"
__all__ = ["Score", "ExtremaResult", "PreparedValues", "compare_scores", "compare_abs_scores", "compare_t", "extrema", "score", "AuditResult", "audit_topk"]
