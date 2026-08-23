.PHONY: sim model test host clean

sim:
	$(MAKE) -C kernel_morbius_plus sim

model:
	$(MAKE) -C host_morbius_plus model

test:
	$(MAKE) -C kernel_morbius_plus sim
	$(MAKE) -C host_morbius_plus test

host:
	$(MAKE) -C host_morbius_plus

clean:
	$(MAKE) -C kernel_morbius_plus clean
	$(MAKE) -C host_morbius_plus clean
