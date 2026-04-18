# coding: utf-8
# @author: R.M.

from constants import m_n

def mu(m_chi, A):
    m_t = A * m_n
    return (m_chi * m_t) / (m_chi + m_t)
    