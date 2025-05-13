ifneq (, $(shell which julia))
.PHONY: dev
dev:
ifneq ($(shell julia --eval 'using MSIAC' 2> /dev/null; echo $$?), 0)
	@echo "Dev MSIAC : add package to Julia dependency"; 
	@julia --eval 'using Pkg; Pkg.develop(path=".")';
endif

.PHONY: test
test: dev
	@echo "Test MSIAC"
	@julia -t $(shell nproc) --eval 'using Pkg; Pkg.test("MSIAC")'

.PHONY: cleanall
cleanall: clean
	@echo "rm MSIAC"
	@julia --eval 'using Pkg; try Pkg.rm("MSIAC") catch throw @info " MSIAC not in Julia distribution. nothing to do here" end'

else
.PHONY: test
test:
	@echo "julia executable not found!"

.PHONY: cleanall
cleanall: clean
endif

.PHONY: clean
clean:
