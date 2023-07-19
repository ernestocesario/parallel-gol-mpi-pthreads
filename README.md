# Parallel Implementation of a cellular automata (Game Of Life) with MPI and pthreads

This is a **parallel implementation** of **John Conway's Game of Life** (done for the **Parallel Algorithms and Distributed Systems** course assignment at university), which uses **MPI** technology to distribute computation across multiple machines (**distributed memory**) and **multithreading** (**shared memory**) to further optimize performance within each individual machine.  
With this combination, the cellular automata can be executed efficiently and quickly on a large scale.  

![alt text](https://github.com/ernestocesario/parallel-gol-mpi-pthreads/blob/main/screenshot.png)

<br>

## General Explanation
### Partitioning
The **domain** of the CA is **partitioned** into submatrices on which each **MPI process** will work.
In turn, each **MPI process** will **partition** its submatrix so that computation can be divided among **multiple threads** within a machine.

### Operation
- Each MPI process partitions its submatrix so that internal threads can work in parallel.
- The threads run a loop until **N_STEP** (which is the value of how many CA steps need to be computed) becomes 0, within which they do:
	- **Sending edges and receiving halo columns** (done only by the main thread of the MPI process)  
	- **Computing each cell** (parallel threads)  
	- **Sending the submatrix to an MPI process that will take care of printing each step** (done only by the main thread of the MPI process)  

## Use
### How to use the Configuration file and the Input file
In the **configuration file** you will find the following **parameters** that you can change:
- **BENCHMARK MODE** (0 - off, 1 - on)
- **N_ROWS**
- **N_COLS**
- **X_PART**
- **Y_PART**
- **N_THREADS**
- **N_STEP**
- **INPUT FILE PATH** (can be empty)

The **input file** is used to give an **initial state** to our **cellular automata**. In this case the file is always assigned to the first MPI worker process.  
You can create an input file (which will be loaded into the array of the first worker process) or use one of those available in the `input_files_ca` folder.  
**Note:** The **input file** must be able to be **contained** in a matrix allocated on an MPI process, which means you must make sure that the **size** (**rows and columns**) of the input file is **less than or equal** to that of a **submatrix** of an **MPI process**.  
**Note:** Separate **columns with a tab** and **rows with a newline**  

### How to Compile and Run
**Install mpich** with: `sudo apt install mpich`  
**Install allegro 4** with: `sudo apt install liballegro4-dev`  
**Compile** the code with: `./compile` *or* `perl ./compile`  
**Run** the code with: `./run N_PROCESS_TO_USE` *or* `perl ./run N_PROCESS_TO_USE`  

<br>

## Detailed explanation
### Y-axis Partitioning (in columns)
**Partitioning** on the **Y-axis** is done by allocating for each **MPI process** a matrix with the **same number of rows** as the domain of the cellular automata (**N_ROWS**) and number of columns equal to the **ratio** of the number of columns of the CA domain (**N_COLS**) to the number of partitions chosen on the Y-axis (**Y_PART**).  
Each MPI process will therefore have a submatrix to compute of size **N_ROWS * (N_COLS / Y_PART)**.  
Actually, the **allocated submatrix** will be **N_ROWS * ((N_COLS / Y_PART) + 2)**.  
Two columns are added because to compute the cells at the edges of our CA we need the state of some cells that are in other processes. These two extra columns then act as our **"halo columns"**, i.e., columns that will receive the state of the CA from the left and right MPI processes (if they exist) to compute the state of the cells at the edges within our MPI process.  
**Note:** It is important that the number of columns of the CA domain (**N_COLS**) is divisible by the number of partitions on the Y axis (**Y_PART**).  

### X-axis Partitioning
Each **MPI process** will **partition** its **submatrix** in order to split the computation between threads, and thus take advantage of parallelism in **shared memory**.  
Partitioning occurs at the beginning in each MPI process.  
All cells in the matrix (including cells in halo columns) are escaped and a simple rule is applied:  
- If it is a cell **belonging to a halo columns** -> **ignore it** (it should not be computed!).  
- If it is a cell **"at the edges"** (i.e. its computation requires consulting a halo cell (cell of a halo columns)) -> **Insert it into a partition of cells at the edges**.  
- If it is an **"inside" cell** (i.e., its computation does NOT require consulting halo cells) -> **Insert it into an inside cell partition**.
  
Doing so creates **X_PART** partitions of inner cells and **X_PART** partitions of outer cells.
The partitions are stored in **two pools** one of **partitions containing only internal cells**, and one of **partitions containing only external cells**.  

### Operation
Each **MPI process** begins by **initializing** all the necessary data structures, and then goes on to start the **threads** that will perform the **CA computation**.  
The **MAXIMUM** number of threads that can be used is equal to the number of partitions on the X-axis (**X_PART**).  
In fact, **having more threads** than **partitions** means leaving the extra threads **idle**, which then only results in a burden, as the system will have to handle them anyway (**unnecessary overhead**). This is why the **platform limits** the maximum number of threads to **X_PART**.
The threads to be started, are equal to N_THREADS - 1, since the main thread also participates in the CA computation.  
  
The **computation** consists of a **series of steps** looped until **N_STEPS** becomes **0**, some of which are done only by one thread:  
- **Sending and receiving halo columns** by NON-blocking operations (**Overlap communications with computations**).  
- **Computation of partitions in the pool of INTERNAL cell partitions** (done in parallel by all)  
- *Wait operation to wait until our process has received halo columns from neighbors*  
- **Computation of partitions in the pool of cell partitions AT BORDER** (done in parallel by all)  
- *Wait operation to wait for "neighbor" processes to have received our halo columns*  
- **Sending the current state of the CA to the MPI process that acts as the display**  

<br>

## Operation of pools
**Pools** have been implemented with the **WorkPool class**.  
They **allow** multiple threads to request partitions of cells to be computed.  
**Increasing** the number of **partitions** will result in **finer granularity** (fewer cells to be computed in each partition), and thus **more parallelism**, but at the same time **increases overhead**, as threads will have to interface many more times with the pool.  
In this case of CA, it is quite useless to have many partitions, since every cell in the cellular automata must be computed in the same way, **but in different cases**, such as computation of a CA that requires the use of **complex functions** for cells that have a certain value, then the use of **dynamic mapping** using a **pool** with **fine granularity increases the degree of parallelism**, compared to a **static mapping** between **partitions** and **threads** (since some threads may have only **"simple"** cells to compute, and thus would end immediately and then **idle** until the threads that have **"complex"** cells have finished their computation).