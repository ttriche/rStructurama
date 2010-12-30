#include "integer.h"
#include <iostream>
#include <climits>
#include <sstream>
#include <cstring>
#include <string>
#include <cmath>

using namespace std;

#if defined(MIN)
#undef MIN
#endif
#if defined(MAX)
#undef MAX
#endif

#define MIN(x,y)  ((x)<(y) ? (x) : (y))
#define MAX(x,y)  ((x)>(y) ? (x) : (y))
#define MAX_COMPONENT  ((unsigned short)(-1))
#define COMPONENT_BITS  (sizeof(unsigned short)*CHAR_BIT)
#define LOG_2_10    3.3219280948873623478703194294894



Integer::Integer(void) {

	createInteger(1);
}

Integer::Integer(string s) {

	createInteger( (int)ceil(LOG_2_10*s.size()/COMPONENT_BITS) );
    Integer digit(1);
    for (int i=0; i<s.size(); i++) 
		{
        multiplySmallInteger(*this, 10, *this);
        digit.c[0] = s[i] - '0';
        addInteger(*this, digit, *this);
		}
}

Integer::Integer(int x) {

	createInteger(x);
}

Integer::~Integer(void) {

	delete [] c;
}

Integer& Integer::operator=(const Integer &a) {

	if ( a.numComponents != numComponents )
		{
		delete [] c;
		createInteger(a.numComponents);
		}
	for (int i=0; i<a.numComponents; i++)
		c[i] = a.c[i];
	return *this;
}

Integer& Integer::operator=(const unsigned long &a) {

	unsigned long val = a;
    unsigned long carry = 0;
	carry += val >> COMPONENT_BITS;
	val &= MAX_COMPONENT;
	
	int numDigits = 1;
	if (carry > 0)
		numDigits = 2;
	
	delete [] c;
	createInteger(numDigits);
	c[0] = val;
	if (numDigits > 1)
		c[1] = carry;
		
	return *this;
}

Integer Integer::operator+(const Integer &a) {

	int n = numComponents;
	int m = a.numComponents;
	Integer result(MAX(n,m)+1);
	addInteger(*this, a, result);
	return result;
}

Integer Integer::operator*(const Integer &a) {

	int n = numComponents;
	int m = a.numComponents;
	Integer result(n+m+1);
	multiplyInteger(*this, a, result);
	return result;
}

Integer Integer::operator*(const unsigned short a) {

	int n = numComponents;
	int m = 1;
	Integer result(n+m+1);
	multiplySmallInteger(*this, a, result);
	return result;
}

ostream &operator<<(ostream &s, Integer &a) {

	string iStr = a.integerToString(a);
	s << iStr;
	return s;
}

void Integer::addInteger(const Integer &left, const Integer &right, Integer &result) {

    unsigned long carry = 0;
    int i;
    for(i=0; i<left.numComponents || i<right.numComponents || carry != 0; i++) 
		{
        unsigned long partial_sum = carry;
        carry = 0;
        if (i < left.numComponents)  
			partial_sum += left.c[i];
        if (i < right.numComponents) 
			partial_sum += right.c[i];
        if (partial_sum > MAX_COMPONENT) 
			{
            partial_sum &= MAX_COMPONENT;
            carry = 1;
			}
        result.c[i] = (unsigned short)partial_sum;
		}
    for ( ; i < result.numComponents; i++) 
		result.c[i] = 0; 
}

void Integer::createInteger(int x) {

	numComponents = x;
	c = new unsigned short[numComponents];
	zeroOutInteger();
}

void Integer::divideSmallInteger(Integer &left, unsigned short right, Integer &result) {

    unsigned long dividend = 0;
    for (int i = left.numComponents - 1; i >= 0; i--) 
		{
        dividend |= left.c[i];
        result.c[i] = dividend/right;
        dividend = (dividend % right) << COMPONENT_BITS;
		}
}

bool Integer::isZeroInteger(Integer &x) {

    for(int i=0; i < x.numComponents; i++) 
		{
        if (x.c[i] != 0) 
			return false;
		}
    return true;
}

double Integer::ln(void) {

	string iStr = integerToString(*this);
	int n = iStr.size();
	iStr.insert(1, ".");
	istringstream buf(iStr);
	double v;
	buf >> v;
	return log(v) + (n-1)*log(10.0);
}

unsigned short Integer::modSmallInteger(Integer &left, unsigned short right) {

    unsigned long modTwoPower = 1;
    unsigned long result = 0;
    for(int i=0; i<left.numComponents; i++) 
		{
        for(int bit=0; bit<COMPONENT_BITS; bit++) 
			{
            if ((left.c[i] & (1 << bit)) != 0) 
				{
                result += modTwoPower;
                if (result >= right) 
                    result -= right;
				}
            modTwoPower <<= 1;
            if (modTwoPower >= right) 
                modTwoPower -= right;
			}
		}
    return (unsigned short)result;
}

void Integer::multiplyInteger(const Integer &left, const Integer &right, Integer &result) {

    unsigned long carry = 0;
    int maxSizeNoCarry;
    int leftMaxComponent  = left.numComponents - 1;
    int rightMaxComponent = right.numComponents - 1;
    maxSizeNoCarry = leftMaxComponent + rightMaxComponent;
	int i;
    for(i=0; i <= maxSizeNoCarry || carry != 0; i++) 
		{
        unsigned long partialSum = carry;
        carry = 0;
        int lidx = MIN(i, leftMaxComponent);
        int ridx = i - lidx;
        while(lidx >= 0 && ridx <= rightMaxComponent) 
			{
            partialSum += ((unsigned long)left.c[lidx])*right.c[ridx];
            carry += partialSum >> COMPONENT_BITS;
            partialSum &= MAX_COMPONENT;
            lidx--; ridx++;
			}
        result.c[i] = partialSum;
		}
    for ( ; i<result.numComponents; i++) 
		result.c[i] = 0; 
}

void Integer::multiplySmallInteger(Integer &left, unsigned short right, Integer &result) {

    unsigned long carry = 0;
    int i;
    for(i=0; i<left.numComponents || carry != 0; i++) 
		{
        unsigned long partial_sum = carry;
        carry = 0;
        if (i < left.numComponents)  
			partial_sum += left.c[i]*right;
        carry = partial_sum >> COMPONENT_BITS;
        result.c[i] = (unsigned short)(partial_sum & MAX_COMPONENT);
		}
    for ( ; i<result.numComponents; i++) 
		result.c[i] = 0; 
}

string Integer::integerToString(Integer &x) {

    int i, result_len;
	char *result = new char[ (int)ceil(COMPONENT_BITS*x.numComponents/LOG_2_10) + 2 ];
    Integer ten(1);
    ten.c[0] = 10;
	Integer val;
	val = x;

    if ( isZeroInteger(val) == true ) 
		{
        strcpy(result, "0");
		} 
	else 
		{
        for (i = 0; !isZeroInteger(val); i++) 
			{
            result[i] = (char)modSmallInteger(val, 10) + '0';
            divideSmallInteger(val, 10, val);
			}
        result[i] = '\0';
		}

    //reverse string
	result_len = strlen(result);
	for(i=0; i < result_len/2; i++) 
		{
		char temp = result[i];
		result[i] = result[result_len - i - 1];
		result[result_len - i - 1] = temp;
		}	
		
	string retStr = result;
	delete [] result;
    return retStr;
}

void Integer::zeroOutInteger(void) {

	for (int i=0; i<numComponents; i++)
		c[i] = 0;
}


