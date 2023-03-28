use crate::util::Average;
use once_cell::sync::Lazy;
use std::collections::HashMap;
use std::sync::Mutex;

pub static FUNCTION_TIMINGS: Lazy<Mutex<HashMap<String, Vec<usize>>>> = Lazy::new(|| {
    let m = HashMap::new();
    Mutex::new(m)
});

pub fn time(function_name: &str, t: usize) {
    let mut locked_local_times = FUNCTION_TIMINGS.lock().unwrap();
    match locked_local_times.get_mut(function_name) {
        Some(times) => times.push(t),
        None => {
            let times: Vec<usize> = vec![t];
            locked_local_times.insert(String::from(function_name), times);
        }
    }
}

pub fn print_timings() {
    let map_guard = FUNCTION_TIMINGS.lock().unwrap();
    let keys = map_guard.keys();
    for key in keys {
        println!("{}: {}", key, map_guard.get(key).unwrap().avg())
    }
}
