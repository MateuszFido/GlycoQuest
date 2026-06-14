use glycoquest::{parse_cli, run};
use std::process;

fn main() {
    let params = match parse_cli(std::env::args()) {
        Ok(params) => params,
        Err(err) => {
            err.print().expect("failed to write error");
            process::exit(err.exit_code());
        }
    };

    process::exit(run(&params));
}
