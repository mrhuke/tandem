# Variable definition
CC = g++
LD = g++
CFLAGS = -g #-DDEBUG
LDFLAGS = -lm 
PROG = tandem
SRCS = tool.cpp gammaTone.cpp filter.cpp feature.cpp voicedMask.cpp tandem.cpp pitch.cpp segmentation.cpp mScaleInten.cpp
OBJS = tool.o gammaTone.o filter.o feature.o voicedMask.o tandem.o pitch.o segmentation.o mScaleInten.o
RM = rm
# end

all:$(PROG) clean

# Link
 $(PROG):$(OBJS)
	$(LD) $(LDFLAGS) $(OBJS) -o $(PROG)

# Compile
%.o:%.cpp
	$(CC) $(CFLAGS) -c $<

clean:
	$(RM) $(OBJS)

