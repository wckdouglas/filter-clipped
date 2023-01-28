use filter_clipped::wrapper;
fn main() {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();
    wrapper();
}
