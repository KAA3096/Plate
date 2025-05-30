import streamlit as st

st.set_page_config(page_title="👋", layout="wide")

st.markdown("""
    <h1 style="text-align:center; font-size: 50px;">Филиал Московского государственного университета имени М. В. Ломоносова в городе Сарове</h1>
""", unsafe_allow_html=True)
st.markdown("""
    <h1 style="text-align:center; font-size: 40px;">Кафедра математики</h1>
""", unsafe_allow_html=True)
st.markdown("""
    <h1 style="text-align:center; font-size: 35px;">Группа ВМ - 124</h1>
""", unsafe_allow_html=True)

# Дополнительное изображение по центру
st.image("logo.jpg", width=300, use_container_width=True)

st.markdown("""
    <h1 style="text-align:center; font-size: 40px;">Численное исследование вихре-потенциального обтекания пластины в канале</h1>
""", unsafe_allow_html=True)

st.markdown("""
    <h1 style="text-align:left; font-size: 35px;">Презентацию подготовили:</h1>
""", unsafe_allow_html=True)

# Данные участников
participants = [
    {"name": "Головня Никита", "photo": "0.jpg", "role": "Построение сеток в Gmsh, генерация геометрии, численное решение (FEniCS) для симметричного обтекания, построение численного решения для симметричного вихре-потенциального обтекания пластины в канале (FEniCS), сопоставление данных функции тока и геометрии в вихревой области с данными по модели Навье-Стокса,
    программная реализация, визуализация, структура презентации."},
    {"name": "Александр Романенко", "photo": "1.jpg", "role": "Верификация численного решения симметричного обтекания при аппроксимации полиномами степени $p = 1,2,3$, при трех значениях циркуляции $\Gamma$,       
при трех значениях полувысоты пластины $h$."},
    {"name": "Гашигуллин Камиль", "photo": "2.jpg", "role": "Верификация численного решения несимметричного обтекания пластины в канале при различных углах наклона, оформление презентации."},
    {"name": "Коврижных Анастасия", "photo": "3.jpg", "role": "Постановка задачи, геометрическая модель, математическая модель, оформление презентации."},
    {"name": "Сержантов Артемий", "photo": "4.jpg", "role": "Переключение слайдов."},
]

# Вывод участников в две строки
row1 = participants[:3]
row2 = participants[3:]

cols1 = st.columns(3)
for i, participant in enumerate(row1):
    with cols1[i]:
        st.image(participant["photo"], width=200)
        st.markdown(f"""
            <h3 style="margin: 0; text-align: left;">{participant['name']}</h3>
            <p style="font-size: 16px; margin: 0; text-align: left;"><i>{participant['role']}</i></p>
        """, unsafe_allow_html=True)

cols2 = st.columns(3)
for i, participant in enumerate(row2):
    with cols2[i]:
        st.image(participant["photo"], width=200)
        st.markdown(f"""
            <h3 style="margin: 0; text-align: left;">{participant['name']}</h3>
            <p style="font-size: 16px; margin: 0; text-align: left;"><i>{participant['role']}</i></p>
        """, unsafe_allow_html=True)

st.markdown("""
    <h2 style="text-align:left;">О презентации</h2>
    <p style="text-align:left; font-size: 18px;">
        Численное исследование вихре-потенциального обтекания пластины в канале представляет 
             значительный интерес для аэродинамики, гидродинамики и смежных инженерных дисциплин. 
             Результаты работы могут быть использованы для оптимизации аэродинамических систем, работающих в ограниченных пространствах,
             таких как элементы гидротурбин, вентиляционные установки и тд.
    </p>
    <p style="text-align:left; font-size: 18px;">
        Данный проект можно посмотреть и скачать на 
        <a href="https://github.com/KAA3096/Plate" target="_blank" style="font-weight: bold;">
        GitHub</a>.
    </p>
""", unsafe_allow_html=True)

col1, col2 = st.columns(2)

with col1:
    st.subheader("🔗 GitHub-репозиторий")
    st.image("qr_github.png", caption="Plate", width=250)
