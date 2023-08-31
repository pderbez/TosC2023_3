/*
 * GField.cpp
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

#include <iostream>
#include "GField.hpp"

GField::GField(unsigned int poly) : m_poly(poly)
{
	if (poly == 2) // Field F_2
	{
		m_dim = 1;
		m_card = 2;
		m_space_mult_table = new GFSymbol [m_card*m_card];
		m_mult_table = new GFSymbol * [m_card];
		m_inv_table = new GFSymbol [m_card];
		for (unsigned int i = 0; i < m_card; ++i) m_mult_table[i] = m_space_mult_table + i*m_card;
		
		m_inv_table[0] = 0;
		m_inv_table[1] = 1;
		
		m_mult_table[0][0] = 0;
		m_mult_table[0][1] = 0;
		m_mult_table[1][0] = 0;
		m_mult_table[1][1] = 1;
	}
	else // Field F_2^m
	{
		m_dim = 0;
		while (poly != 1)
		{
			poly >>= 1;
			++m_dim;
		}
		m_card = 1 << m_dim;
		
		m_space_mult_table = new GFSymbol [m_card*m_card];
		m_mult_table = new GFSymbol * [m_card];
		m_inv_table = new GFSymbol [m_card];
		
		for (unsigned int i = 0; i < m_card; ++i) m_mult_table[i] = m_space_mult_table + i*m_card;
		
		for (unsigned int i = 0; i < m_card; ++i)
		{
			for (unsigned int j = i; j < m_card; ++j)
			{
				unsigned int x = i;
				unsigned int y = j;
				
				// Compute and store i*j
				GFSymbol res = 0;
				for (unsigned int k = 0; k < m_dim; ++k)
				{
					if (x%2 == 1) res ^= y;
					x >>= 1;
					y <<= 1;
					if (y >= m_card) y ^= m_poly;
				}
				m_mult_table[i][j] = res;
				if (j != i) m_mult_table[j][i] = res;
			}
		}
		
		// Construction of the table of inverses  
		m_inv_table[0] = 0;
		unsigned int cpt_inv = 0;
		for (unsigned int i = 1; i < m_card; ++i)
		{
			for (unsigned int j = i; j < m_card; ++j)
			{
				if (m_mult_table[i][j] == 1)
				{
					m_inv_table[i] = j;
					++cpt_inv;
					if (i!=j)
					{
						m_inv_table[j] = i;
						++cpt_inv;
					}
					break;
				}
			}
		}
		
		// Check if each non-zero element has an inverse
		if (cpt_inv != m_card - 1)
		{
			std::cout << "Error in Field construction -> it's not a field!!!" << std::endl;
			exit(EXIT_FAILURE);
		}
	}
}


GField::GField(GField const& gf) : m_poly(gf.m_poly)
{
	m_dim = gf.m_dim;
	m_card = gf.m_card;
	
	m_space_mult_table = new GFSymbol [m_card*m_card];
	m_mult_table = new GFSymbol * [m_card];
	m_inv_table = new GFSymbol [m_card];
	
	{
		const unsigned int bound = m_card*m_card;
		for (unsigned int i = 0; i < bound; ++i) m_space_mult_table[i] = gf.m_space_mult_table[i];
	}
	
	{
		const unsigned int bound = m_card;
		for (unsigned int i = 0; i < bound; ++i) m_mult_table[i] = m_space_mult_table + i*bound;
	}
	
	{
		const unsigned int bound = m_card;
		for (unsigned int i = 0; i < bound; ++i) m_inv_table[i] = gf.m_inv_table[i];
	}
}

GField::~GField()
{
	delete[] m_space_mult_table;
	delete[] m_mult_table;
	delete[] m_inv_table;
}

