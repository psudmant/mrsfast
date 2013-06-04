mrsFAST 2.5.0.4 clone with pipe support
---------------------------------------

This version of mrsFAST is a clone of the [mrsFAST v2.5.0.4 ](http://mrsfast.sourceforge.net/Home) mapping software by Faraz Hach.

In the version the Reads.c code has been patched. Instead of traversing the file to count the lines and then seeking to 0 to read in again, this program reads the file into a buffer and then reads the buffer back into memory. 

A simpler set of read functions is also implemented.

These changes allow mrsfast to read from pipes, very helpful.

FUTURE IMPROVEMENTS
------------------

Instead of reading into a buffer: contruct the file on the fly - allocate memory for say, 500,000 reads, then realloc as necessary up, or down.

contact psudmant@gmail.com for questions/comments
