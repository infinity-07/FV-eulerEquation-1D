#./bin/program
 
# 编译器
CC = mpicxx

# 编译选项
CFLAGS = -Wall -O3 

# 目录路径
src = ./src
inc = ./include
build = ./build
bin = ./bin

# 源文件
SRCS = $(src)/main.cpp $(src)/WENOFV.cpp

# 头文件路径
INCLUDES = -I$(inc)

# 目标文件 (将 .cpp 转为 .o，但输出到 build 文件夹)
OBJS = $(SRCS:$(src)/%.cpp=$(build)/%.o)

# 可执行文件目标
TARGET = $(bin)/program

all: $(TARGET)

# 生成最终可执行文件
$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^

# 生成 .o 文件，输出到 build 文件夹
$(build)/%.o: $(src)/%.cpp
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@


clean:
	rm -f $(build)/*.o $(TARGET)
	