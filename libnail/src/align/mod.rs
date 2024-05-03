pub mod structs;

mod cloud_search;
pub use cloud_search::{
    cloud_search_backward, cloud_search_forward, prune_and_scrub, scrub_co_located,
    CloudSearchParams, CloudSearchScores,
};

mod forward;
pub use forward::forward;

mod backward;
pub use backward::backward;

mod posterior;
pub use posterior::posterior;

mod optimal_accuracy;
pub use optimal_accuracy::optimal_accuracy;

mod scoring;
pub use scoring::{e_value, null_one_score, null_two_score, p_value, Bits, Nats, Score};

mod traceback;
pub use traceback::traceback;

mod needleman_wunsch;
pub use needleman_wunsch::{needleman_wunsch, SimpleTraceStep};
