This is my own implementation of std::deque. 
The main methods of the specified containers are implemented, as well as support for iterators and reverse iterators.
The deque is implemented as an array of pointers to blocks of elements, where the blocks can be located far apart in 
memory.
All methods are strong exception safety. Additionally, it supports move semantics, allowing for efficient transfer of 
objects between containers without the need for copying. Deque can also be configured using allocators, enabling 
control over the allocation and deallocation of memory for the container's elements.
