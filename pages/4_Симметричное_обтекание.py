import streamlit as st
import base64
import pandas as pd

menu = st.sidebar.radio(
    '***',
    (
        "Аппроксимация полиномами степени $p=1,2,3$",
        "Сравнение результатов",
        "Циркулляция",
        "Полувысота",
        "Программная реализация",
    )
)

def load_gif(gif_path):
    try:
        with open(gif_path, "rb") as file_:
            contents = file_.read()
            data_url = base64.b64encode(contents).decode("utf-8")
            return st.markdown(
                f'<img src="data:image/gif;base64,{data_url}" width="100%">',
                unsafe_allow_html=True,
            )
    except FileNotFoundError:
        st.error(f"Файл {gif_path} не найден!")
    except Exception as e:
        st.error(f"Ошибка загрузки GIF: {str(e)}")

if menu == "Аппроксимация полиномами степени $p=1,2,3$":
    r"""
    ##### Численное решение задачи вихре-потенциального симметричного обтекания пластины в канале при аппроксимации полиномами степени $p=1,2,3$ 
    """
    
    mesh_size = st.selectbox(
        "Размер сетки",  
        ["985", "3738", "14788"]  
    )
    
    poly_degree = st.selectbox(
        "Степень полинома",  
        ["p = 1", "p = 2", "p = 3"]  
    )

    if mesh_size == "985":
        if poly_degree == "p = 1":
            load_gif("/home/anastasia/Документы/rep/985_p1.gif")
            r"""
            -----------------------------------------------
            * Число ячеек сетки: 985
            * Число узлов сетки: 549
            * Число искомых дискретных значений: 549
            * Макс. завихренность: 2.801862
            * Мин. завихренность: -19.281812
            * Итерации: 5
            * $\Gamma = -1$
            * $\omega_0 = -2.0$
            ------------------------------------------------
            """
            
        elif poly_degree == "p = 2":
            load_gif("/home/anastasia/Документы/rep/985_p2.gif")
            r"""
            -----------------------------------------------
            * Число ячеек сетки: 985
            * Число узлов сетки: 549
            * Число искомых дискретных значений: 2082
            * Макс. завихренность: 5.189690
            * Мин. завихренность: -21.665007
            * Итерации: 6
            * $\Gamma = -1$
            * $\omega_0 = -2.0$
            ------------------------------------------------
            """
            
        elif poly_degree == "p = 3":
            load_gif("/home/anastasia/Документы/rep/985_p3.gif")
            r"""
            -----------------------------------------------
            * Число ячеек сетки: 985
            * Число узлов сетки: 549
            * Число искомых дискретных значений: 4600
            * Макс. завихренность: 5.114822
            * Мин. завихренность: -17.977754
            * Итерации: 6
            * $\Gamma = -1$
            * $\omega_0 = -2.0$
            ------------------------------------------------
            """
    if mesh_size == "3738":
        if poly_degree == "p = 1":
            load_gif("/home/anastasia/Документы/rep/3738_p1.gif")
            r"""
            -----------------------------------------------
            * Число ячеек сетки: 3738
            * Число узлов сетки: 1980
            * Число искомых дискретных значений: 1980
            * Макс. завихренность: 3.111834
            * Мин. завихренность: -17.913237
            * Итерации: 5
            * $\Gamma = -1$
            * $\omega_0 = -2.0$
            ------------------------------------------------
            """
            
        elif poly_degree == "p = 2":
            load_gif("/home/anastasia/Документы/rep/3738_p2.gif")
            r"""
            -----------------------------------------------
            * Число ячеек сетки: 3738
            * Число узлов сетки: 1980
            * Число искомых дискретных значений: 17152
            * Макс. завихренность: 5.632535
            * Мин. завихренность: -20.340337
            * Итерации: 16
            * $\Gamma = -1$
            * $\omega_0 = -2.0$
            ------------------------------------------------
            """
            
        elif poly_degree == "p = 3":
            load_gif("/home/anastasia/Документы/rep/3738_p3.gif")
            r"""
            -----------------------------------------------
            * Число ячеек сетки: 3738
            * Число узлов сетки: 1980
            * Число искомых дискретных значений: 17152
            * Макс. завихренность: 5.632822
            * Мин. завихренность: -20.340337
            * Итерации: 16
            * $\Gamma = -1$
            * $\omega_0 = -2.0$
            ------------------------------------------------
            """
            
    if mesh_size == "14788":
        if poly_degree == "p = 1":
            load_gif("/home/anastasia/Документы/rep/14788_p1.gif")
            r"""
            -----------------------------------------------
            * Число ячеек сетки: 14788
            * Число узлов сетки: 7615
            * Число искомых дискретных значений: 7615
            * Макс. завихренность: 3.133065
            * Мин. завихренность: -16.742390
            * Итерации: 16
            * $\Gamma = -1$
            * $\omega_0 = -2.0$
            ------------------------------------------------
            """
            
        elif poly_degree == "p = 2":
            load_gif("/home/anastasia/Документы/rep/14788_p2.gif")
            r"""
            -----------------------------------------------
            * Число ячеек сетки: 14788
            * Число узлов сетки: 7615
            * Число искомых дискретных значений: 30017
            * Макс. завихренность: 6.126937
            * Мин. завихренность: -18.957333
            * Итерации: 12
            * $\Gamma = -1$
            * $\omega_0 = -2.0$
            ------------------------------------------------
            """
            
        elif poly_degree == "p = 3":
            load_gif("/home/anastasia/Документы/rep/14788_p3.gif")
            r"""
            -----------------------------------------------
            * Число ячеек сетки: 14788
            * Число узлов сетки: 7615
            * Число искомых дискретных значений: 67207
            * Макс. завихренность: 5.340547
            * Мин. завихренность: -18.714660
            * Итерации: 36
            * $\Gamma = -1$
            * $\omega_0 = -2.0$
            ------------------------------------------------
            """
if menu == "Сравнение результатов":
    r"""
    ##### Число степеней свободы в зависимости от размера сетки и порядка аппроксимации  $p$
    """
    data = {
        "Размер сетки": [985, 3738, 14788],
        "Ст. свободы": [549, 1980, 7615],
        "Макс. завихр.": [2.80, 3.11, 3.13],
        "Мин. завихр.": [-19.28, -17.91, -16.74]
    }

    df_p1 = pd.DataFrame({
        "Размер сетки": data["Размер сетки"],
        "Степени свободы": data["Ст. свободы"],
        "Макс. завихренность": data["Макс. завихр."],
        "Мин. завихренность": data["Мин. завихр."]
    })
    
    data = {    
        "Размер сетки": [985, 3738, 14788],
        "Ст. свободы": [2082, 17152, 30017],
        "Макс. завихр.": [5.19, 5.63, 6.13],
        "Мин. завихр.": [-21.67, -20.34, -18.96]
    }

    df_p2 = pd.DataFrame({
        "Размер сетки": data["Размер сетки"],
        "Степени свободы": data["Ст. свободы"],
        "Макс. завихренность": data["Макс. завихр."],
        "Мин. завихренность": data["Мин. завихр."]
    })
    
    data = { 
        "Размер сетки": [985, 3738, 14788],
        "Ст. свободы": [4600, 17152, 67207],
        "Макс. завихр.": [5.11, 5.63, 5.34],
        "Мин. завихр.": [-17.98, -20.34, -18.71]
    }

    df_p3 = pd.DataFrame({
        "Размер сетки": data["Размер сетки"],
        "Степени свободы": data["Ст. свободы"],
        "Макс. завихренность": data["Макс. завихр."],
        "Мин. завихренность": data["Мин. завихр."]
    })


    r"""Порядок аппроксимации $p = 1$"""
    st.dataframe(
        df_p1,
        hide_index=True,
        use_container_width=True,
        column_config={
            "Размер сетки": st.column_config.NumberColumn(width="small"),
            "Степени свободы": st.column_config.NumberColumn(width="medium"),
            "Макс. завихренность": st.column_config.NumberColumn(format="%.2f"),
            "Мин. завихренность": st.column_config.NumberColumn(format="%.2f")
        }
    )

    r"""Порядок аппроксимации $p = 2$"""
    st.dataframe(
        df_p2,
        hide_index=True,
        use_container_width=True,
        column_config={
            "Размер сетки": st.column_config.NumberColumn(width="small"),
            "Степени свободы": st.column_config.NumberColumn(width="medium"),
            "Макс. завихренность": st.column_config.NumberColumn(format="%.2f"),
            "Мин. завихренность": st.column_config.NumberColumn(format="%.2f")
        }
    )

    r"""Порядок аппроксимации $p = 3$"""
    st.dataframe(
        df_p3,
        hide_index=True,
        use_container_width=True,
        column_config={
            "Размер сетки": st.column_config.NumberColumn(width="small"),
            "Степени свободы": st.column_config.NumberColumn(width="medium"),
            "Макс. завихренность": st.column_config.NumberColumn(format="%.2f"),
            "Мин. завихренность": st.column_config.NumberColumn(format="%.2f")
        }
    )
        
            
if menu == "Циркулляция":
    r"""
    ##### Численное решение задачи вихре-потенциального симметричного обтекания пластины в канале при различной циркуляции $\Gamma$ 
    """
    gamma = st.selectbox(
        "Значение циркуляции",  
        ["-1", "-2", "-3", "-4", "-5", "-6", "-7", "-8", "-9"]  
    )
    
    if gamma == "-1":
        load_gif("/home/anastasia/Документы/rep/gamma_1.gif")
        r"""
            -----------------------------------------------
            * Число ячеек сетки: 366098
            * Число узлов сетки: 184150
            * Число искомых дискретных значений: 184150
            * Макс. завихренность: 3.211237
            * Мин. завихренность: -16.231253
            * Итераций: 17
            * $\omega_0 =$ -1
            * p = 1
            ------------------------------------------------
        """
    if gamma == "-2":
        load_gif("/home/anastasia/Документы/rep/gamma_2.gif")
        r"""
            -----------------------------------------------
            * Число ячеек сетки: 366098
            * Число узлов сетки: 184150
            * Число искомых дискретных значений: 184150
            * Макс. завихренность: 3.005262
            * Мин. завихренность: -15.305808
            * Итераций: 280
            * $\omega_0 =$ -2
            * p = 1
            ------------------------------------------------
        """
    if gamma == "-3":
        load_gif("/home/anastasia/Документы/rep/gamma_3.gif")
        r"""
            -----------------------------------------------
            * Число ячеек сетки: 366098
            * Число узлов сетки: 184150
            * Число искомых дискретных значений: 184150
            * Макс. завихренность: 2.812799
            * Мин. завихренность: -14.185900
            * Итераций: 147
            * $\omega_0 =$ -3
            * p = 1
            ------------------------------------------------
        """
    if gamma == "-4":
        load_gif("/home/anastasia/Документы/rep/gamma_4.gif")
        r"""
            -----------------------------------------------
            * Число ячеек сетки: 366098
            * Число узлов сетки: 184150
            * Число искомых дискретных значений: 184150
            * Макс. завихренность: 2.636768
            * Мин. завихренность: -13.782692
            * Итераций: 109
            * $\omega_0 =$ -4
            * p = 1
            ------------------------------------------------
        """
    if gamma == "-5":
        load_gif("/home/anastasia/Документы/rep/gamma_5.gif")
        r"""
            -----------------------------------------------
            * Число ячеек сетки: 366098
            * Число узлов сетки: 184150
            * Число искомых дискретных значений: 184150
            * Макс. завихренность: 2.511733
            * Мин. завихренность: -12.690229
            * Итераций: 86
            * $\omega_0 =$ -5
            * p = 1
            ------------------------------------------------
        """
    if gamma == "-6":
        load_gif("/home/anastasia/Документы/rep/gamma_6.gif")
        r"""
            -----------------------------------------------
            * Число ячеек сетки: 366098
            * Число узлов сетки: 184150
            * Число искомых дискретных значений: 184150
            * Макс. завихренность: 2.424755
            * Мин. завихренность: -12.354141
            * Итераций: 73
            * $\omega_0$ = -6
            * p = 1
            ------------------------------------------------
        """
    if gamma == "-7":
        load_gif("/home/anastasia/Документы/rep/gamma_7.gif")
        r"""
            -----------------------------------------------
            * Число ячеек сетки: 366098
            * Число узлов сетки: 184150
            * Число искомых дискретных значений: 184150
            * Макс. завихренность: 2.364655
            * Мин. завихренность: -11.941278
            * Итераций: 69
            * $\omega_0 =$ -7
            * p=1
        """
    if gamma == "-8":
        load_gif("/home/anastasia/Документы/rep/gamma_8.gif")
        r"""
            -----------------------------------------------
            * Число ячеек сетки: 366098
            * Число узлов сетки: 184150
            * Число искомых дискретных значений: 184150
            * Макс. завихренность: 2.331714
            * Мин. завихренность: -11.825124
            * Итераций: 17
            * $\omega_0 =$ -8
            * p=1
            ------------------------------------------------
        """
    if gamma == "-9":
        load_gif("/home/anastasia/Документы/rep/gamma_9.gif")
        r"""
            -----------------------------------------------
            * Число ячеек сетки: 366098
            * Число узлов сетки: 184150
            * Число искомых дискретных значений: 184150
            * Макс. завихренность: 2.312357
            * Мин. завихренность: -11.721477
            * Итераций: 65
            * $\omega_0 =$ -9
            * p = 1
            ------------------------------------------------
        """
        
if menu == "Полувысота":
    r"""
    ##### Численное решение задачи вихре-потенциального симметричного обтекания пластины в канале при изменении полувысоты пластины $h$
    """
    
    height = st.selectbox(
        "Полувысота пластины",  
        ["0.1", "0.25", "0.6"]  
    )
    
    if height == "0.1":
        load_gif("/home/anastasia/Документы/rep/height_01.gif")
        r"""
            -----------------------------------------------
            * Число ячеек сетки: 366098
            * Число узлов сетки: 184150
            * Число искомых дискретных значений: 184150
            * Макс. завихренность: 8.913647
            * Мин. завихренность: -46.390642
            * Итераций: 23
            * $\Gamma =$ -1.5
            * $\omega_0 =$ -1.5
            * p = 1
            ------------------------------------------------
        """
    if height == "0.25":
        load_gif("/home/anastasia/Документы/rep/height_025.gif")
        r"""
            -----------------------------------------------
            * Число ячеек сетки: 366098
            * Число узлов сетки: 184150
            * Число искомых дискретных значений: 184150
            * Макс. завихренность: 4.079369
            * Мин. завихренность: -20.610520
            * Итераций: 70
            * $\Gamma =$ -2.0
            * $\omega_0 =$ -2.0
            * p = 1
            ------------------------------------------------
        """
    if height == "0.6":
        load_gif("/home/anastasia/Документы/rep/height_06.gif")
        r"""
            -----------------------------------------------
            * Число ячеек сетки: 366098
            * Число узлов сетки: 184150
            * Число искомых дискретных значений: 184150
            * Макс. завихренность: 2.322699
            * Мин. завихренность: -11.729287
            * Итераций: 86
            * $\Gamma =$ -11.0
            * $\omega_0 =$ -11.0
            * p = 1
            ------------------------------------------------
        """
    
if menu == "Программная реализация":
    r"""
    ##### Программная реализация
    """
    
    with st.expander("Определение функционального пространства"):
        code = """
            mesh = Mesh("plate.xml")
            boundaries = MeshFunction("size_t", mesh, "plate_facet_region.xml")
            ds = Measure("ds", subdomain_data=boundaries)

            deg = 1
            V = FunctionSpace(mesh, "CG", deg)
        """
        st.code(code, language="python")
        
    with st.expander("Граничные условия"):
        code = """
        u_l = Expression("x[1]", degree=deg)
        u_r = Expression("x[1]", degree=deg)
        u_t = Expression("x[1]", degree=deg)
        u_b = Expression("0.0", degree=deg)

        bcs1 = [DirichletBC(V, u_l, boundaries, 2),
                DirichletBC(V, u_r, boundaries, 3),
                DirichletBC(V, u_t, boundaries, 4),
                DirichletBC(V, u_b, boundaries, 5)]

        """
        st.code(code, language="python")
        
    with st.expander("Начальные условия"):
        code = """
        Gamma = -2.0
        omega_0 = '-2.0'
        omega_k = Function(V)
        omega_k.interpolate(Expression(omega_0, degree=deg))
        u_k = Function(V)

        """
        st.code(code, language="python")
        
    with st.expander("Итерационный процесс"):
        code = """
        max_iter = 500
        tolerance = 1e-12

        for k in range(max_iter):
            print(f"\n== Итерация {k} ==")
            a = dot(grad(u), grad(TestFunction(V))) * dx
            L = omega_k * TestFunction(V) * dx
            solve(a == L, u_k, bcs1)


            indicator = conditional(lt(u_k, 0.0), 1.0, 0.0)
            area = assemble(indicator * dx)

            if near(area, 0.0):
                print("Область ψ < 0 исчезла, остановка итераций")
        """
        st.code(code, language="python")
    with st.expander("Проверка сходимости"):
        code = """
       if k > 0:
        change = errornorm(u_k, u_prev, 'L2')
        print(f"Change: {change}")
        if change < tolerance:
            print("Достигнута сходимость")
            break
        """
        st.code(code, language="python")
           
    
        


        
    
    
        


    


            
 
            
            
            
            
            
            
            
