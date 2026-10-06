"""ExtremaRank: exact studentized subset extrema and shared-donor rank audit."""
from .core import Score, ExtremaResult, PreparedValues, compare_scores, compare_abs_scores, compare_t, extrema, score
from .audit import AuditResult, audit_topk
from .unpaired import WelchScore, WelchAuditResult, PreparedWelch, compare_welch, audit_welch
from .diagnostics import PreparedPaired
from .robustness import audit_feature_robustness

__version__ = "0.4.0"
__all__ = ["Score", "ExtremaResult", "PreparedValues", "compare_scores", "compare_abs_scores", "compare_t", "extrema", "score", "AuditResult", "audit_topk", "WelchScore", "WelchAuditResult", "PreparedWelch", "compare_welch", "audit_welch", "PreparedPaired", "audit_feature_robustness"]
