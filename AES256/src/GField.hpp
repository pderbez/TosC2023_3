/*
 * GField.h
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


#ifndef DEF_GFIELD
#define DEF_GFIELD

#include <cstdint>

typedef uint8_t GFSymbol;

/**
 *  @brief     Represent a finite field GF(2^n)
 *  @warning   The polynomial has to be irreductible!
 */

class GField
{
	public:
		GField(unsigned int poly);
		GField(GField const& gf);
		~GField();

		unsigned int getDim() const {return m_dim;};
		unsigned int getCard() const {return m_card;};
		unsigned int getPoly() const {return m_poly;};

		GFSymbol inverse(const GFSymbol val) const {return m_inv_table[val];};
		GFSymbol multiply(const GFSymbol val1, const GFSymbol val2) const {return m_mult_table[val1][val2];};

	private:

		const unsigned int m_poly;
		unsigned int m_dim;
		unsigned int m_card;

		GFSymbol * m_space_mult_table;
		GFSymbol ** m_mult_table;
		GFSymbol * m_inv_table;
};


#endif
