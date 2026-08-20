.PHONY: sim model test host clean

sim:
	$(MAKE) -C hw sim

model:
	$(MAKE) -C sw model

test:
	$(MAKE) -C hw sim
	$(MAKE) -C sw test

host:
	$(MAKE) -C sw

clean:
	$(MAKE) -C hw clean
	$(MAKE) -C sw clean
