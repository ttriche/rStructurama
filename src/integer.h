#ifndef INTEGER_H
#define INTEGER_H

#include <string>
#include <iostream>

using namespace std;



class Integer {

		   friend ostream   &operator<<(ostream &s, Integer &a);
		   
	public:
	                        Integer(void);
                            Integer(string s);
							Integer(int x);
							~Integer(void);
				  Integer   &operator=(const Integer &a);
				  Integer   &operator=(const unsigned long &a);
				  Integer   operator+(const Integer &a);
				  Integer   operator*(const Integer &a);
				  Integer   operator*(const unsigned short a);
				   double   ln(void);

	private:
	                 void   addInteger(const Integer &left, const Integer &right, Integer &result);
	                 void   createInteger(int x);
					 void   divideSmallInteger(Integer &left, unsigned short right, Integer &result);
				   string   integerToString(Integer &x);
					 bool   isZeroInteger(Integer &x);
		   unsigned short   modSmallInteger(Integer &left, unsigned short right);
		             void   multiplyInteger(const Integer &left, const Integer &right, Integer &result);
					 void   multiplySmallInteger(Integer &left, unsigned short right, Integer &result);
					 void   zeroOutInteger(void);
		   unsigned short   *c;
					  int   numComponents;
};

#endif



