PYTHON := python2.7
PIP := ~/.local/bin/pip2.7
VIRTUALENV := ~/.local/bin/virtualenv
SANDBOX := .python2.7

.DEFAULT_GOAL := requeirements

$(PIP):
	TEMP=`mktemp /tmp/get-pip.XXXXXX.py`;\
	curl -o $$TEMP https://bootstrap.pypa.io/get-pip.py;\
	$(PYTHON) $$TEMP --user;\
	rm $$TEMP

$(VIRTUALENV): $(PIP)
	$(PIP) install virtualenv --user

$(SANDBOX): $(VIRTUALENV)
	$(VIRTUALENV) -p $(PYTHON) $@

.PHONY: requeirements
requeirements: $(SANDBOX)
	source $(SANDBOX)/bin/activate;\
	pip install -r requirements.txt
