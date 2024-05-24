use std::sync::atomic::{AtomicUsize, Ordering};

static ID_COUNTER: AtomicUsize = AtomicUsize::new(0);

#[derive(PartialEq, Eq, Hash, Clone, Copy, Debug)]
pub struct Id(usize);

/// Produce an `Id` using the global Id counter.
///
/// This is thread-safe, as the counter is maintained using an `AtomicUsize`.
pub fn next_id() -> Id {
    let id_num = ID_COUNTER.fetch_add(1, Ordering::SeqCst);
    Id(id_num)
}

/// Produce the last Id that was generated with `id::next_id()`.
pub fn last_id() -> Id {
    let id_num = ID_COUNTER.load(Ordering::SeqCst);
    Id(id_num)
}
