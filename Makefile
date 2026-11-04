# Releases are published by .github/workflows/publish.yml when a GitHub release is made,
# so this file has no upload target.

.PHONY: help
help: ## This help message
	@echo -e "$$(grep -hE '^\S+:.*##' $(MAKEFILE_LIST) | sed -e 's/:.*##\s*/:/' -e 's/^\(.\+\):\(.*\)/\\x1b[36m\1\\x1b[m:\2/' | column -c2 -t -s :)"

.PHONY: install
install: ## Install gwr-inversion with every extra (editable), with the dev group
	pip install -e ".[all]" --group dev --config-settings editable_mode=compat

.PHONY: lint
lint: ## Run ruff and mypy the way CI does
	ruff check gwr_inversion/
	mypy gwr_inversion/

.PHONY: test
test: lint ## Run lint, then pytest
	python -m pytest test/ -q

.PHONY: build
build: clean ## Build the wheel and sdist into an empty dist/
	python -m build

.PHONY: clean
clean: ## Remove build artifacts
	rm -rf build/ dist/ gwr_inversion.egg-info/
