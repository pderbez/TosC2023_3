/*
 * GFElement.h
 * 
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 * 
 * 
 */


#ifndef DEF_GFELEMENT
#define DEF_GFELEMENT

#include <iostream>
#include <iomanip>
#include "GField.hpp"

/**
 *  @brief     Represent an element of a finite field GF(2^n)
 *  @warning   It is assumed that operators take elements that belong to the same field.
 */

class GFElement
{
	public:
		GFElement() {};
		GFElement(const GFSymbol value) : m_value(value) {};
		GFElement(GFElement const & element) : m_value(element.m_value) {};
		
		GFElement& operator=(const GFElement element);
		GFElement& operator=(const GFSymbol value);
		
		~GFElement() {};
		
		void inverse();
		GFElement getInverse() const;
#if (GF_POLY != 0x02)
		static GFSymbol inverse(const GFSymbol value) {return m_gf.inverse(value);};
#else
		static GFSymbol inverse(const GFSymbol value) {return value;};
#endif
		static unsigned int getDim() { return m_gf.getDim(); };
		static unsigned int getCard() { return m_gf.getCard(); };
		
		GFElement& operator+=(const GFElement element);
		GFElement& operator-=(const GFElement element);
		GFElement& operator*=(const GFElement element);
		GFElement& operator/=(const GFElement element);
		
		GFElement& operator+=(const GFSymbol value);
		GFElement& operator-=(const GFSymbol value);
		GFElement& operator*=(const GFSymbol value);
		GFElement& operator/=(const GFSymbol value);
		
	
	private:
		GFSymbol m_value;
		
		const static GField m_gf; 
	
		
	friend std::ostream& operator<<( std::ostream &flux, const GFElement var);
	
	friend GFElement operator+(const GFElement element1, const GFElement element2);
	friend GFElement operator-(const GFElement element1, const GFElement element2);
	friend GFElement operator*(const GFElement element1, const GFElement element2);
	friend GFElement operator/(const GFElement element1, const GFElement element2);
	
	friend GFElement operator+(const GFElement element, const GFSymbol val);
	friend GFElement operator-(const GFElement element, const GFSymbol val);
	friend GFElement operator*(const GFElement element, const GFSymbol val);
	friend GFElement operator/(const GFElement element, const GFSymbol val);
	
	friend GFElement operator+(const GFSymbol val, const GFElement element);
	friend GFElement operator-(const GFSymbol val, const GFElement element);
	friend GFElement operator*(const GFSymbol val, const GFElement element);
	friend GFElement operator/(const GFSymbol val, const GFElement element);
	
	friend bool operator==(const GFElement element1, const GFElement element2) { return (element1.m_value == element2.m_value); }; 	
	friend bool operator==(const GFElement element, const GFSymbol val) { return (element.m_value == val); }; //useful to try == 1 or == 0
	friend bool operator==(const GFSymbol val, const GFElement element) { return (element.m_value == val); }; 
	
	friend bool operator!=(const GFElement element1, const GFElement element2) { return (element1.m_value != element2.m_value); }; 	
	friend bool operator!=(const GFElement element, const GFSymbol val) { return (element.m_value != val); }; 
	friend bool operator!=(const GFSymbol val, const GFElement element) { return (element.m_value != val); };
	
	friend bool operator<(const GFElement element1, const GFElement element2) { return (element1.m_value < element2.m_value); }; // to use std::map, std::set
};



// Operators

inline GFElement& GFElement::operator=(const GFElement element) 
{
    m_value = element.m_value;
    return *this; 
}

inline GFElement& GFElement::operator=(const GFSymbol value)
{
	m_value = value;
	return *this;
}

inline GFElement& GFElement::operator+=(const GFElement element)
{
	m_value ^= element.m_value;
	return *this;
}

inline GFElement& GFElement::operator-=(const GFElement element)
{
	m_value ^= element.m_value;
	return *this;
}

inline GFElement& GFElement::operator+=(const GFSymbol value)
{
	m_value ^= value;
	return *this;
}

inline GFElement& GFElement::operator-=(const GFSymbol value)
{
	m_value ^= value;
	return *this;
}

#if (GF_POLY != 0x02) // not the finite field F_2

inline GFElement& GFElement::operator*=(const GFElement element)
{
	m_value = m_gf.multiply(m_value,element.m_value);
	return *this;
}

inline GFElement& GFElement::operator/=(const GFElement element)
{
	m_value = m_gf.multiply(m_value,m_gf.inverse(element.m_value));
	return *this;
}

inline GFElement& GFElement::operator*=(const GFSymbol value)
{
	m_value = m_gf.multiply(value,m_value);
	return *this;
}

inline GFElement& GFElement::operator/=(const GFSymbol value)
{
	m_value = m_gf.multiply(m_gf.inverse(value),m_value);
	return *this;
}

#else // the finite field F_2

inline GFElement& GFElement::operator*=(const GFElement element)
{
	m_value &= element.m_value;
	return *this;
}

inline GFElement& GFElement::operator/=(const GFElement element)
{
	m_value &= element.m_value;
	return *this;
}

inline GFElement& GFElement::operator*=(const GFSymbol value)
{
	m_value &= value;
	return *this;
}

inline GFElement& GFElement::operator/=(const GFSymbol value)
{
	m_value &= value;
	return *this;
}

#endif



// Member Functions 

#if (GF_POLY != 0x02) // not the finite field F_2

inline GFElement GFElement::getInverse() const
{
	return GFElement(m_gf.inverse(m_value));
}

inline void GFElement::inverse()
{
	m_value = m_gf.inverse(m_value);
}

#else

inline GFElement GFElement::getInverse() const
{
	return *this;
}

inline void GFElement::inverse()
{
}

#endif


// Friend Functions 

inline GFElement operator+(const GFElement element1, const GFElement element2)
{
	return GFElement(element1.m_value ^ element2.m_value);
}

inline GFElement operator-(const GFElement element1, const GFElement element2)
{
	return GFElement(element1.m_value ^ element2.m_value);
}

inline GFElement operator+(const GFElement element, const GFSymbol val)
{
	return GFElement(element.m_value ^ val);
}

inline GFElement operator-(const GFElement element, const GFSymbol val)
{
	return GFElement(element.m_value ^ val);
}
	
inline GFElement operator+(const GFSymbol val, const GFElement element)
{
	return GFElement(element.m_value ^ val);
}

inline GFElement operator-(const GFSymbol val, const GFElement element)
{
	return GFElement(element.m_value ^ val);
}

#if (GF_POLY != 0x02) // not the finite field F_2

inline GFElement operator*(const GFElement element1, const GFElement element2)
{
	return GFElement(GFElement::m_gf.multiply(element1.m_value,element2.m_value));
}

inline GFElement operator/(const GFElement element1, const GFElement element2)
{
	return GFElement(GFElement::m_gf.multiply(GFElement::m_gf.inverse(element2.m_value), element1.m_value));
}

inline GFElement operator*(const GFElement element, const GFSymbol val)
{
	return GFElement(GFElement::m_gf.multiply(val,element.m_value));
}

inline GFElement operator/(const GFElement element, const GFSymbol val)
{
	return GFElement(GFElement::m_gf.multiply(GFElement::m_gf.inverse(val),element.m_value));
}

inline GFElement operator*(const GFSymbol val, const GFElement element)
{
	return GFElement(GFElement::m_gf.multiply(val,element.m_value));
}

inline GFElement operator/(const GFSymbol val, const GFElement element)
{
	return GFElement(GFElement::m_gf.multiply(val,GFElement::m_gf.inverse(element.m_value)));
}


#else

inline GFElement operator*(const GFElement element1, const GFElement element2)
{
	return GFElement(element1.m_value & element2.m_value);
}

inline GFElement operator/(const GFElement element1, const GFElement element2)
{
	return GFElement(element1.m_value & element2.m_value);
}

inline GFElement operator*(const GFElement element, const GFSymbol val)
{
	return GFElement(element.m_value & val);
}

inline GFElement operator/(const GFElement element, const GFSymbol val)
{
	return GFElement(element.m_value & val);
}

inline GFElement operator*(const GFSymbol val, const GFElement element)
{
	return GFElement(element.m_value & val);
}

inline GFElement operator/(const GFSymbol val, const GFElement element)
{
	return GFElement(element.m_value & val);
}

#endif

inline std::ostream& operator<<( std::ostream &flux, const GFElement var)
{
	flux << std::hex << std::setfill('0') << std::setw(2) << static_cast<unsigned int>(var.m_value) << std::dec;
	return flux;
}




#endif
