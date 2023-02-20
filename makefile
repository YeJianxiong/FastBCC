GCC=g++ -O3 -fcilkplus
OBJS=FastBCC.o
FastBCC: $(OBJS)
	$(GCC) -o FastBCC $(OBJS)
clean:
	rm $(OBJS); rm FastBCC
FastBCC.o: FastBCC.cpp CycleTimer.h
	$(GCC) -c FastBCC.cpp -lpthread -g -DDISCRETE 
