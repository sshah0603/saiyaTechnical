# Grading entry points. Three targets, exactly as specified in the brief:
#
#   make setup      install dependencies
#   make pipeline   database, Part 1, then all tables and plots for Parts 2-4
#   make dashboard  start the local dashboard server
#
# Everything runs with the interpreter that invoked make, so this works the
# same inside a virtualenv, a Codespace, or a bare container.

PYTHON ?= python3
PORT   ?= 8501

.DEFAULT_GOAL := help
.PHONY: help setup pipeline dashboard test report clean

help:
	@echo "make setup      Install dependencies from requirements.txt"
	@echo "make pipeline   Build the database and generate all outputs"
	@echo "make dashboard  Start the dashboard on port $(PORT)"
	@echo "make test       Run the test suite"
	@echo "make clean      Remove the database and generated outputs"

# Plain install first. Some images mark the system Python as externally
# managed (PEP 668), which rejects a plain install; the fallback handles that
# without forcing a virtualenv on the grader.
setup:
	-$(PYTHON) -m pip install --upgrade pip
	$(PYTHON) -m pip install -r requirements.txt \
	  || $(PYTHON) -m pip install --break-system-packages -r requirements.txt

# Part 1 initializes the schema and loads the CSV. run_pipeline.py then
# produces every table and plot for Parts 2, 3 and 4 into outputs/.
pipeline:
	$(PYTHON) load_data.py
	$(PYTHON) run_pipeline.py

# python -m streamlit rather than the bare command, so this works whether or
# not pip's script directory is on PATH.
dashboard:
	$(PYTHON) -m streamlit run app.py --server.port $(PORT) --server.address 0.0.0.0

test:
	$(PYTHON) -m pytest -q

report:
	$(PYTHON) make_report.py

clean:
	rm -f cell_counts.db report.html
	rm -rf outputs __pycache__ .pytest_cache
