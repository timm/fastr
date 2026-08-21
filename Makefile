include .dot/Makefile

$(shell mkdir -p ~/tmp)

Cpu=10
What=contrast.py
Flag=-halve
Data=$(HOME)/gits/moot/optimize

xargs: ## run What+Flag on every Data csv, Cpu at a time
	@find $(Data) -name '*.csv' | \
	xargs -P $(Cpu) -n 1 -I{} python3 $(What) --file {} $(Flag)

~/tmp/contrast.txt: ## save the sorted xargs run here
	@$(MAKE) -s xargs | tee $@
