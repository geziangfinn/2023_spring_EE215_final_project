core_router: core_router.cpp# Build the core_router
	g++ -o core_router ../src/core_router.cpp

run: # Run the core router for all benchmarks
	make run_phase_1
	make run_phase_2

run_phase_1: # Run the core router for phase_1 benchmarks
	./core_router ../benchmarks bench1 ../out
	./core_router ../benchmarks bench2 ../out
	./core_router ../benchmarks bench3 ../out
	./core_router ../benchmarks bench4 ../out
	
run_phase_2: # Run the core router for phase_2 benchmarks
	./core_router ../benchmarks bench5 ../out_preferred_routing
	./core_router ../benchmarks fract2 ../out_preferred_routing
	./core_router ../benchmarks industry1 ../out_preferred_routing
	./core_router ../benchmarks primary1 ../out_preferred_routing

.PHONY:
clean: # Clean the build directory, and routing results
	rm -f core_router
	rm -f ../out/*.route