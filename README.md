//////////
This version of mrsfast has a patched Reads.c file. 
Instead of traversing the file to count the lines and then seeking 
to 0 to read in again, this program reads the file into a buffer 
and then reads the buffer back into memory. 

A simpler set of read functions is also implemented.

These changes allow mrsfast to read from pipes, very helpful.

FUTURE IMPROVEMENTS
instead of reading into a buffer, 
contruct the file on the fly - allocate memory for say, 
500,000 reads, then realloc as necessary up, or down.

-----------------
contact psudmant [ a t ] gmail . com for details
----------------
