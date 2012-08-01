#include <iostream>
#include "bit_vector.h"

bit_vector::bit_vector(int64_t num_bits)
    : _num_bits(num_bits), _db(NULL)
{
  uint64_t num_bytes = _num_bits / 8;
  if (_num_bits % 8 != 0)
  {
    num_bytes += 1;
  }

  if ((size_t)num_bytes < num_bytes)
  {
    std::cout << "Error: won't be able to calloc " << num_bytes << " bytes on architecture with sizeof(size_t) " << sizeof(size_t) << "\n";
  }

  if ((_db = calloc(num_bytes, 1)) == NULL)
  {
    perror("calloc bit_vector()");
    std::cout << "num_bytes " << num_bytes << "\n";
    exit(1);
  }
}

bit_vector::~bit_vector()
{
  free(_db);
}

