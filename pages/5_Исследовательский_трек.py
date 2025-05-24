import streamlit as st
import base64

menu = st.sidebar.radio(
    '***',
    (
        "Модель Навье-Стокса",
        "Сравнение по функции тока в вихревой области",
        "Сравнение по геометрии вихревой зоны",
        "Программная реализация"
    )
)
if menu == "Модель Навье-Стокса":
    r"""
    ##### Модель Навье-Стокса 
    
    **Уравнения Навье-Стокса** для несжимаемой жидкости:
 

    $\begin{aligned}\frac{\partial \mathbf{u}}{\partial t} + (\mathbf{u} \cdot \nabla) \mathbf{u} = 
    -\frac{1}{\rho} \nabla p + \nu \nabla^2 \mathbf{u}\end{aligned}$ 
  
    
     $\begin{aligned}\nabla \cdot \mathbf{u} = 0
    \end{aligned}$
    
    где:
    
    $\mathbf{u} - $вектор скорости;
    
    $p - $ давление;
    
    $\rho - $ плотность жидкости;
    
    $\nu - $ кинематическая вязкость.
    
    **Масштабирование:**
    
    Для обобщения результатов введем безразмерные переменные:
    
    - пространственные координаты: $\bar{x} = x / H$, $\bar{y} = y / H$
    - скорость: $\bar{\mathbf{u}} = \mathbf{u} / U$
    - время: $\bar{t} = H^2 / \nu$
    - давление: $\bar{p} = p H / (\mu U)$
    
    где $U -$ характерная скорость на входе, $H -$ пространственный масштаб;
    
    После масштабирования уравнения принимают вид:
    
     $\begin{aligned}\frac{\partial \bar{\mathbf{u}}}{\partial \bar{t}} + {\operatorname{Re}}(\bar{\mathbf{u}} \cdot \nabla) \bar{\mathbf{u}} = 
    -\nabla \bar{p} + \nabla^2 \bar{\mathbf{u}}\end{aligned}$
    
    $\begin{aligned}\nabla \cdot \bar{\mathbf{u}} = 0\end{aligned}$
    
    где $\operatorname{Re} = \frac{\rho U L}{\mu}$ — число Рейнольдса, $\mu -$ динамическая вязкость, $\nu = \mu / \rho$.
    
    """
if menu == "Сравнение по функции тока в вихревой области":
    r"""
    ##### Сравнение расчетных данных функции тока в вихревой области с данными по модели Навье-Стокса 
    """
if menu == "Сравнение по геометрии вихревой зоны":
    r"""
    ##### Сравнение расчетных данных по геометрии вихревой зоны с данными по модели Навье-Стокса 
    """
if menu == "Программная реализация":
    r"""
    ##### Программная реализация 
    """
    with st.expander("Определение функциональных пространств"):
        code = """
            # Пространства
            V = VectorFunctionSpace(mesh, 'P', 1)  # Пространство для скорости
            Q = FunctionSpace(mesh, 'P', 1)         # Пространство для давления
        """
        st.code(code, language="python")
        
    with st.expander("Граничные условия"):
        code = """
            # Граничные условия для скорости
            inflow_profile = Expression(('1.0', '0.0'), degree=2)
            bcu_inflow = DirichletBC(V, inflow_profile, boundaries, 2)  # Вход
            bcu_walls  = DirichletBC(V, Constant((1, 0)), boundaries, 4)  # Стенки
            bcu_plate  = DirichletBC(V, Constant((0, 0)), boundaries, 5)  # Пластина
            bcu = [bcu_inflow, bcu_walls, bcu_plate]

            # Граничное условие для давления (на выходе)
            bcp_outflow = DirichletBC(Q, Constant(0), boundaries, 3)
            bcp = [bcp_outflow]
        """
        st.code(code, language="python")
        
    with st.expander("Определение форм и функций"):
        code = """
            u = TrialFunction(V)  # Пробная функция для скорости
            v = TestFunction(V)   # Тестовая функция для скорости
            p = TrialFunction(Q)  # Пробная функция для давления
            q = TestFunction(Q)   # Тестовая функция для давления

            u_n = Function(V)     # Скорость на предыдущем шаге
            u_  = Function(V)     # Текущая скорость
            p_n = Function(Q)     # Давление на предыдущем шаге
            p_  = Function(Q)     # Текущее давление
        """
        st.code(code, language="python")
        
    with st.expander("Формулировка слабых форм уравнений"):
        code = """
            # Тензор напряжений
            def sigma(u, p):
                return 2*mu*epsilon(u) - p*Identity(len(u))

            # Слабая форма для уравнения Навье-Стокса
            F1 = rho*dot((u - u_n)/k, v)*dx + ... 
            a1 = lhs(F1)
            L1 = rhs(F1)

            # Слабая форма для уравнения Пуассона (давление)
            a2 = dot(nabla_grad(p), nabla_grad(q))*dx
            L2 = ... 

            # Слабая форма для коррекции скорости
            a3 = dot(u, v)*dx
            L3 = ... 
        """
        st.code(code, language="python")
        
    with st.expander("Сборка матриц и применение граничных условий"):
        code = """
            A1 = assemble(a1)
            A2 = assemble(a2)
            A3 = assemble(a3)
            [bc.apply(A1) for bc in bcu]  # Применение граничных условий для скорости
            [bc.apply(A2) for bc in bcp]  # Применение граничных условий для давления
        """
        st.code(code, language="python")
        
    with st.expander("Временной цикл"):
        code = """
            for n in range(num_steps):
                t += dt
                # Решение для скорости
                solve(A1, u_.vector(), b1, 'bicgstab', 'hypre_amg')
                # Решение для давления
                solve(A2, p_.vector(), b2, 'bicgstab', 'hypre_amg')
                # Коррекция скорости
                solve(A3, u_.vector(), b3, 'cg', 'sor')
                # Сохранение результатов
                xdmffile_u.write(u_, t)
                timeseries_u.store(u_.vector(), t)
                # Обновление переменных
                u_n.assign(u_)
                p_n.assign(p_)
        """
        st.code(code, language="python")
        
        
