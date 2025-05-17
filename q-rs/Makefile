SHELL := /bin/bash

test: fmt lint
	cargo test

fmt:
	cargo fmt

lint:
	cargo clippy

rustup:
	curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh

update:
	rustup update
	cargo update

check:
	cargo check

publish:
	cargo login
	cargo publish --dry-run