#ifndef __BIT_VECTOR_HH
#define __BIT_VECTOR_HH

#include <iostream>

using namespace std;

class bit_vector
{
public:
  bit_vector(int64_t _num_bits);
  ~bit_vector();

  void set_one(int64_t);

  void set_zero(int64_t);
  
  bool get(int64_t) const;

  int64_t num_bits_set() const;
  
  void print() const;

protected:
  int64_t _num_bits;
  void *_db;
};

inline
void bit_vector::print() const
{
	for (int64_t i = _num_bits-1; i >= 0; i--) {
		cout << get(i);
		if (i % 8 == 0)
			cout << " ";
	}
	cout << endl;
}

inline
int64_t bit_vector::num_bits_set() const
{
	return _num_bits;
}

inline
void bit_vector::set_one(int64_t bit_idx)
{
  static_cast<unsigned char*>(_db)[bit_idx >> 3] |= (1 << (bit_idx & 7));
}

inline
void bit_vector::set_zero(int64_t bit_idx)
{
  static_cast<unsigned char*>(_db)[bit_idx >> 3] &= ~(1 << (bit_idx & 7));
}
  
inline
bool bit_vector::get(int64_t bit_idx) const
{
  return (static_cast<unsigned char*>(_db)[bit_idx >> 3] >> (bit_idx & 7)) & 1;
}

#endif

