PYTHON := python2.6
PIP := ~/.local/bin/pip2.6
VIRTUALENV := ~/.local/bin/virtualenv
SANDBOX := .sandbox

.DEFAULT_GOAL := requeirements

$(PIP):
	TEMP=`mktemp /tmp/get-pip.XXXXXX.py`;\
	curl -o $$TEMP https://bootstrap.pypa.io/get-pip.py;\
	$(PYTHON) $$TEMP --user;\
	rm $$TEMP

$(VIRTUALENV): $(PIP)
	$(PIP) install virtualenv --user

$(SANDBOX): $(VIRTUALENV)
	$(VIRTUALENV) -p $(PYTHON) --system-site-packages $@

.PHONY: requeirements
requeirements: $(SANDBOX)
	source $(SANDBOX)/bin/activate;\
	pip install -r requirements.txt
