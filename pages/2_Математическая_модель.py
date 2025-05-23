import streamlit as st
from PIL import Image

r"""
    ##### Математическая модель
    
    **Вихрь скорости**
        
    $\begin{aligned}
    \nabla \times \bm{v}= \{0, 0, \omega\},
    \quad \omega = \frac{\partial v_2}{\partial x_1} - \frac{\partial v_1}{\partial x_2}
    \end{aligned}$   
         
    **Уравнение для функции тока**
    
    Стационарные течения идеальной жидкости
 
    $\begin{aligned}
    {-} \nabla^2 \psi = \omega(\psi),
    \quad \bm x \in \Omega
    \end{aligned}$   
    
    - $\omega = 0 - $ потенциальное течение
    - $\omega \neq 0 - $ вихревое течение
     
    **Краевые условия**
    
    Равномерный поток (скорость $\bm u = \{1, 0\}$) на входе в канал

    $\begin{aligned}
    \psi(\bm x) = H x_2,
    \quad \bm x \in \Gamma_3
    \end{aligned}$
    
    На выходе (плоскопараллельное течение)
    
    $\begin{aligned}
    \frac{\partial \psi}{\partial n}(\bm x)  = 0, 
    \quad \bm x \in  \Gamma_4
    \end{aligned}$
    
    На границе пластины и линии симметрии
  
    $\begin{aligned}
    \psi(\bm x) = 0,
    \quad \bm x \in \Gamma_1
    \end{aligned}$  
    
    На стенках канала $\psi = \operatorname{const} - $ условие непротекания
    
    $\begin{aligned}
    \psi(\bm x) = H ,
    \quad \bm x \in \Gamma_2
    \end{aligned}$    
    
    **Склейка потенциальных и вихревых течений**
    
    Области течений
    
    - $\psi > 0 - $ потенциальное течение
    - $\psi < 0 - $ вихревое течение
    
    Вихревые течения
    
    - $\omega(\psi) = \omega = \operatorname{const} - $ модель постоянной завихренности
    - $\omega < 0 - $ вихревое течение
    
    Всегда есть есть решение $-$ чисто потенциальное течение
    
    Выделение вихре-потенциального решения по заданной циркуляции
    
    $\begin{aligned}
    \int_{\psi < 0} \omega(\psi) d \bm x = \Gamma
    \end{aligned}$    

    """
