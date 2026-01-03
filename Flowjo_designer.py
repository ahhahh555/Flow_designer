# app.py
import streamlit as st
import pandas as pd
import numpy as np
from typing import Dict, List, Optional
import json
from dataclasses import dataclass, field, asdict
from enum import Enum
import base64
from datetime import datetime
import io
import altair as alt  # 使用Altair代替Plotly

# 设置页面配置
st.set_page_config(
    page_title="流式染色矩阵设计器",
    page_icon="🔬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# 自定义CSS
st.markdown("""
<style>
    .main-header {
        font-size: 2.5rem;
        color: #1E3A8A;
        text-align: center;
        margin-bottom: 2rem;
        font-weight: bold;
    }
    .section-header {
        font-size: 1.8rem;
        color: #2563EB;
        margin-top: 2rem;
        margin-bottom: 1rem;
        border-bottom: 3px solid #3B82F6;
        padding-bottom: 0.5rem;
    }
    .stButton > button {
        background-color: #3B82F6;
        color: white;
        font-weight: bold;
        border-radius: 8px;
        padding: 0.5rem 1rem;
        border: none;
    }
    .stButton > button:hover {
        background-color: #1D4ED8;
        color: white;
    }
    .success-box {
        background-color: #D1FAE5;
        border: 1px solid #10B981;
        border-radius: 8px;
        padding: 1rem;
        margin: 1rem 0;
    }
    .info-box {
        background-color: #DBEAFE;
        border: 1px solid #3B82F6;
        border-radius: 8px;
        padding: 1rem;
        margin: 1rem 0;
    }
    .warning-box {
        background-color: #FEF3C7;
        border: 1px solid #F59E0B;
        border-radius: 8px;
        padding: 1rem;
        margin: 1rem 0;
    }
    .antibody-card {
        background: linear-gradient(135deg, #667eea20 0%, #764ba240 100%);
        border-left: 5px solid #667eea;
        padding: 1rem;
        border-radius: 8px;
        margin: 0.5rem 0;
    }
    .tube-card {
        background: linear-gradient(135deg, #f093fb20 0%, #f5576c40 100%);
        border-left: 5px solid #f5576c;
        padding: 1rem;
        border-radius: 8px;
        margin: 0.5rem 0;
    }
    .control-tube {
        background: linear-gradient(135deg, #4ade8020 0%, #22c55e40 100%);
        border-left: 5px solid #22c55e;
    }
    .experiment-tube {
        background: linear-gradient(135deg, #3b82f620 0%, #2563eb40 100%);
        border-left: 5px solid #2563eb;
    }
</style>
""", unsafe_allow_html=True)

class AntibodyType(Enum):
    """抗体类型"""
    SURFACE = "表面抗体"
    INTRACELLULAR = "胞内抗体"
    VIABILITY = "死活染料"
    FC_BLOCK = "Fc阻断剂"
    OTHER = "其他"

@dataclass
class Antibody:
    """抗体信息"""
    name: str
    short_name: str = ""
    fluorochrome: str = ""
    target: str = ""
    clone: str = ""
    concentration: float = 0.0  # μg/mL
    recommended_use: float = 0.0  # μg/10⁶ cells
    type: AntibodyType = AntibodyType.SURFACE
    catalog_number: str = ""
    lot_number: str = ""
    storage: str = "4°C避光"
    notes: str = ""
    
    def __post_init__(self):
        if not self.short_name:
            # 自动生成简称
            words = self.name.split()
            if words:
                # 取最后一个单词或第一个单词的前几个字母
                self.short_name = words[-1][:8]
            else:
                self.short_name = self.name[:8]

@dataclass
class TubeConfiguration:
    """管子配置"""
    name: str
    description: str
    antibodies: List[str] = field(default_factory=list)
    needs_fixation: bool = False
    is_control: bool = False
    control_type: str = ""  # FMO, Isotype, Single, Blank

def init_session():
    """初始化session state"""
    defaults = {
        'antibodies': {},
        'tubes': {},
        'current_project': f"Flow_Project_{datetime.now().strftime('%Y%m%d_%H%M')}",
        'volumes': {
            'per_tube': 100.0,
            'intracellular_per_tube': 50.0,
            'cell_count': 1.0,  # 10⁶ cells
            'extra_tubes': 2
        }
    }
    
    for key, value in defaults.items():
        if key not in st.session_state:
            st.session_state[key] = value

def load_standard_antibodies():
    """加载标准抗体"""
    standard_antibodies = [
        Antibody(
            name="TruStain FcX™ (anti-mouse CD16/32) Antibody",
            short_name="FcX",
            fluorochrome="None",
            target="CD16/32",
            clone="93",
            concentration=500.0,
            recommended_use=1.0,
            type=AntibodyType.FC_BLOCK,
            catalog_number="101320",
            lot_number="B123456",
            notes="Fc受体阻断剂"
        ),
        Antibody(
            name="Live/Dye eF780",
            short_name="LiveDye",
            fluorochrome="eF780",
            target="Viability",
            clone="N/A",
            concentration=1000.0,
            recommended_use=0.5,
            type=AntibodyType.VIABILITY,
            catalog_number="65-0865-14",
            lot_number="123456",
            notes="死活染料，建议1:1000稀释"
        ),
        Antibody(
            name="BB515 Rat Anti-Mouse CD45",
            short_name="CD45",
            fluorochrome="BB515",
            target="CD45",
            clone="30-F11",
            concentration=200.0,
            recommended_use=0.25,
            type=AntibodyType.SURFACE,
            catalog_number="564590",
            lot_number="789012",
            notes="白细胞标记"
        ),
        Antibody(
            name="α-SMA AF647",
            short_name="α-SMA",
            fluorochrome="AF647",
            target="α-SMA",
            clone="1A4",
            concentration=200.0,
            recommended_use=0.5,
            type=AntibodyType.INTRACELLULAR,
            catalog_number="561847",
            lot_number="345678",
            notes="胞内染色，需固定破膜"
        )
    ]
    
    for ab in standard_antibodies:
        st.session_state.antibodies[ab.name] = ab

def load_standard_tubes():
    """加载标准管子配置"""
    standard_tubes = {
        "Blank": TubeConfiguration(
            name="Blank",
            description="未染色对照，调节电压",
            is_control=True,
            control_type="Blank"
        ),
        "FcX_Only": TubeConfiguration(
            name="FcX_Only",
            description="仅Fc阻断对照",
            antibodies=["TruStain FcX™ (anti-mouse CD16/32) Antibody"],
            is_control=True,
            control_type="Single"
        ),
        "Live_Only": TubeConfiguration(
            name="Live_Only",
            description="死活染料单阳（补偿）",
            antibodies=["TruStain FcX™ (anti-mouse CD16/32) Antibody", "Live/Dye eF780"],
            is_control=True,
            control_type="Single"
        ),
        "CD45_Only": TubeConfiguration(
            name="CD45_Only",
            description="CD45单阳（补偿）",
            antibodies=["TruStain FcX™ (anti-mouse CD16/32) Antibody", "BB515 Rat Anti-Mouse CD45"],
            is_control=True,
            control_type="Single"
        ),
        "αSMA_Only": TubeConfiguration(
            name="αSMA_Only",
            description="α-SMA单阳（补偿，需破膜）",
            antibodies=["TruStain FcX™ (anti-mouse CD16/32) Antibody", "α-SMA AF647"],
            needs_fixation=True,
            is_control=True,
            control_type="Single"
        ),
        "FMO_αSMA": TubeConfiguration(
            name="FMO_αSMA",
            description="荧光减一对照（用于α-SMA设门）",
            antibodies=["TruStain FcX™ (anti-mouse CD16/32) Antibody", "Live/Dye eF780", "BB515 Rat Anti-Mouse CD45"],
            is_control=True,
            control_type="FMO"
        ),
        "Full_Stain": TubeConfiguration(
            name="Full_Stain",
            description="全染实验管",
            antibodies=["TruStain FcX™ (anti-mouse CD16/32) Antibody", "Live/Dye eF780", "BB515 Rat Anti-Mouse CD45", "α-SMA AF647"],
            needs_fixation=True,
            is_control=False
        )
    }
    
    st.session_state.tubes = standard_tubes

def display_antibody_card(antibody: Antibody):
    """显示抗体卡片"""
    type_colors = {
        "表面抗体": "#3B82F6",
        "胞内抗体": "#10B981",
        "死活染料": "#F59E0B",
        "Fc阻断剂": "#EF4444",
        "其他": "#6B7280"
    }
    
    color = type_colors.get(antibody.type.value, "#6B7280")
    
    with st.container():
        col1, col2 = st.columns([4, 1])
        with col1:
            st.markdown(f"**{antibody.name}**")
            st.caption(f"靶标: {antibody.target} | 荧光: {antibody.fluorochrome} | 克隆: {antibody.clone}")
        with col2:
            st.markdown(f'<span style="color:{color}; font-weight:bold;">{antibody.type.value}</span>', unsafe_allow_html=True)
        
        st.markdown(f"浓度: {antibody.concentration} μg/mL | 用量: {antibody.recommended_use} μg/10⁶ cells")
        if antibody.notes:
            st.info(antibody.notes)
        
        st.divider()

def display_tube_card(tube: TubeConfiguration):
    """显示管子卡片"""
    with st.container():
        card_class = "control-tube" if tube.is_control else "experiment-tube"
        
        col1, col2, col3 = st.columns([3, 2, 1])
        with col1:
            st.markdown(f"**{tube.name}**")
            st.caption(tube.description)
        with col2:
            if tube.is_control:
                st.markdown(f'**对照类型:** {tube.control_type}')
            if tube.needs_fixation:
                st.markdown("**需固定破膜**")
        with col3:
            st.markdown(f"**抗体数:** {len(tube.antibodies)}")
        
        if tube.antibodies:
            with st.expander("查看抗体列表"):
                for ab_name in tube.antibodies:
                    if ab_name in st.session_state.antibodies:
                        ab = st.session_state.antibodies[ab_name]
                        st.markdown(f"- {ab.short_name} ({ab.fluorochrome})")
        
        st.divider()

def render_dashboard():
    """仪表盘页面"""
    st.markdown('<div class="main-header">🔬 流式染色矩阵设计器</div>', unsafe_allow_html=True)
    
    # 项目概览卡片
    col1, col2, col3, col4 = st.columns(4)
    
    with col1:
        st.metric("抗体数量", len(st.session_state.antibodies))
    with col2:
        st.metric("管子配置", len(st.session_state.tubes))
    with col3:
        surface_count = sum(1 for ab in st.session_state.antibodies.values() 
                          if ab.type == AntibodyType.SURFACE)
        st.metric("表面抗体", surface_count)
    with col4:
        intracel_count = sum(1 for ab in st.session_state.antibodies.values() 
                           if ab.type == AntibodyType.INTRACELLULAR)
        st.metric("胞内抗体", intracel_count)
    
    # 快速开始指南
    st.markdown("### 🚀 快速开始")
    
    guide_col1, guide_col2 = st.columns(2)
    
    with guide_col1:
        st.info("#### 新手模式")
        if st.button("📥 加载标准实验模板", use_container_width=True):
            load_standard_antibodies()
            load_standard_tubes()
            st.success("✅ 标准模板加载成功！")
            st.rerun()
    
    with guide_col2:
        st.info("#### 自定义模式")
        if st.button("🆕 开始自定义设计", use_container_width=True):
            st.session_state.page = "antibodies"
            st.rerun()
    
    # 功能卡片
    st.markdown("### 🛠️ 功能模块")
    
    func_col1, func_col2, func_col3 = st.columns(3)
    
    with func_col1:
        if st.button("🧪 抗体管理", use_container_width=True):
            st.session_state.page = "antibodies"
            st.rerun()
    
    with func_col2:
        if st.button("🧫 管子设计", use_container_width=True):
            st.session_state.page = "tubes"
            st.rerun()
    
    with func_col3:
        if st.button("🔢 染色矩阵", use_container_width=True):
            st.session_state.page = "matrix"
            st.rerun()
    
    # 最近活动
    if st.session_state.antibodies:
        st.markdown("### 📊 当前项目概览")
        
        overview_col1, overview_col2 = st.columns(2)
        
        with overview_col1:
            st.markdown("#### 最近添加的抗体")
            for ab in list(st.session_state.antibodies.values())[-3:]:
                st.markdown(f"- **{ab.short_name}**: {ab.target} ({ab.fluorochrome})")
        
        with overview_col2:
            if st.session_state.tubes:
                st.markdown("#### 最近添加的管子")
                for tube in list(st.session_state.tubes.values())[-3:]:
                    st.markdown(f"- **{tube.name}**: {len(tube.antibodies)}种抗体")

def render_antibody_management():
    """抗体管理页面"""
    st.markdown('<div class="section-header">🧪 抗体库管理</div>', unsafe_allow_html=True)
    
    tab1, tab2 = st.tabs(["📝 添加抗体", "📋 抗体列表"])
    
    with tab1:
        st.markdown("### 添加新抗体")
        
        with st.form("add_antibody_form"):
            col1, col2 = st.columns(2)
            
            with col1:
                name = st.text_input("抗体全名*", placeholder="例如: TruStain FcX™ Antibody")
                target = st.text_input("靶标", placeholder="例如: CD45")
                fluorochrome = st.text_input("荧光染料", placeholder="例如: AF647")
                clone = st.text_input("克隆号", placeholder="例如: 30-F11")
            
            with col2:
                concentration = st.number_input("浓度 (μg/mL)*", min_value=0.0, value=200.0)
                recommended_use = st.number_input("用量 (μg/10⁶ cells)*", min_value=0.0, value=0.5)
                antibody_type = st.selectbox(
                    "抗体类型*",
                    options=[t.value for t in AntibodyType],
                    index=0
                )
            
            catalog_number = st.text_input("货号")
            lot_number = st.text_input("批号")
            notes = st.text_area("备注")
            
            submitted = st.form_submit_button("✅ 添加抗体", use_container_width=True)
            
            if submitted:
                if not name:
                    st.error("抗体名称不能为空！")
                else:
                    # 转换类型字符串为枚举
                    type_map = {t.value: t for t in AntibodyType}
                    
                    antibody = Antibody(
                        name=name,
                        target=target,
                        fluorochrome=fluorochrome,
                        clone=clone,
                        concentration=concentration,
                        recommended_use=recommended_use,
                        type=type_map[antibody_type],
                        catalog_number=catalog_number,
                        lot_number=lot_number,
                        notes=notes
                    )
                    
                    st.session_state.antibodies[name] = antibody
                    st.success(f"✅ 已成功添加抗体: {name}")
    
    with tab2:
        st.markdown(f"### 抗体库 ({len(st.session_state.antibodies)}种)")
        
        if not st.session_state.antibodies:
            st.info("暂无抗体数据")
        else:
            search_term = st.text_input("🔍 搜索抗体", placeholder="输入名称、靶标或荧光染料搜索")
            
            filtered_antibodies = []
            for ab in st.session_state.antibodies.values():
                if (not search_term or 
                    search_term.lower() in ab.name.lower() or 
                    search_term.lower() in ab.target.lower() or 
                    search_term.lower() in ab.fluorochrome.lower()):
                    filtered_antibodies.append(ab)
            
            if filtered_antibodies:
                for ab in filtered_antibodies:
                    col1, col2 = st.columns([4, 1])
                    with col1:
                        display_antibody_card(ab)
                    with col2:
                        if st.button("删除", key=f"del_{ab.name}", type="secondary"):
                            del st.session_state.antibodies[ab.name]
                            st.rerun()
            else:
                st.warning("未找到匹配的抗体")

def render_tube_design():
    """管子设计页面"""
    st.markdown('<div class="section-header">🧫 实验管子设计</div>', unsafe_allow_html=True)
    
    tab1, tab2 = st.tabs(["🎯 创建管子", "📊 管子列表"])
    
    with tab1:
        st.markdown("### 创建新管子配置")
        
        with st.form("add_tube_form"):
            col1, col2 = st.columns(2)
            
            with col1:
                tube_name = st.text_input("管子名称*", placeholder="例如: Full_Stain")
                description = st.text_area("描述", placeholder="例如: 全染实验管")
            
            with col2:
                col_type1, col_type2 = st.columns(2)
                with col_type1:
                    is_control = st.checkbox("是否为对照管")
                with col_type2:
                    needs_fixation = st.checkbox("需要固定破膜")
                
                if is_control:
                    control_type = st.selectbox(
                        "对照类型",
                        options=["FMO", "Isotype", "Single", "Blank"],
                        index=0
                    )
                else:
                    control_type = ""
            
            # 选择抗体
            if st.session_state.antibodies:
                st.markdown("### 选择抗体")
                selected_antibodies = []
                
                # 按类型分组显示抗体
                antibody_types = {}
                for ab in st.session_state.antibodies.values():
                    if ab.type.value not in antibody_types:
                        antibody_types[ab.type.value] = []
                    antibody_types[ab.type.value].append(ab)
                
                for type_name, antibodies in antibody_types.items():
                    with st.expander(f"{type_name} ({len(antibodies)}种)"):
                        for ab in antibodies:
                            if st.checkbox(f"{ab.name} ({ab.fluorochrome})", key=f"tube_{tube_name}_{ab.name}"):
                                selected_antibodies.append(ab.name)
            else:
                st.warning("请先添加抗体")
                selected_antibodies = []
            
            submitted = st.form_submit_button("✅ 创建管子配置", use_container_width=True)
            
            if submitted:
                if not tube_name:
                    st.error("管子名称不能为空！")
                else:
                    tube = TubeConfiguration(
                        name=tube_name,
                        description=description,
                        antibodies=selected_antibodies,
                        needs_fixation=needs_fixation,
                        is_control=is_control,
                        control_type=control_type
                    )
                    
                    st.session_state.tubes[tube_name] = tube
                    st.success(f"✅ 已成功创建管子: {tube_name}")
    
    with tab2:
        st.markdown(f"### 管子配置 ({len(st.session_state.tubes)}种)")
        
        if not st.session_state.tubes:
            st.info("暂无管子配置")
        else:
            for tube in st.session_state.tubes.values():
                col1, col2 = st.columns([4, 1])
                with col1:
                    display_tube_card(tube)
                with col2:
                    if st.button("删除", key=f"del_tube_{tube.name}", type="secondary"):
                        del st.session_state.tubes[tube.name]
                        st.rerun()

def render_matrix():
    """染色矩阵页面"""
    st.markdown('<div class="section-header">🔢 染色矩阵</div>', unsafe_allow_html=True)
    
    if not st.session_state.tubes or not st.session_state.antibodies:
        st.warning("请先配置抗体和管子")
        return
    
    # 创建矩阵数据
    matrix_data = []
    tube_names = list(st.session_state.tubes.keys())
    antibody_names = list(st.session_state.antibodies.keys())
    
    for tube_name in tube_names:
        tube = st.session_state.tubes[tube_name]
        row = {"管子名称": tube_name, "描述": tube.description}
        
        for ab_name in antibody_names:
            row[ab_name] = "✓" if ab_name in tube.antibodies else "○"
        
        row["固定破膜"] = "✓" if tube.needs_fixation else ""
        row["对照类型"] = tube.control_type if tube.is_control else "实验管"
        matrix_data.append(row)
    
    # 创建DataFrame
    df = pd.DataFrame(matrix_data)
    
    # 显示矩阵
    st.markdown("### 染色矩阵表")
    
    # 使用st.dataframe显示可交互表格
    st.dataframe(
        df,
        use_container_width=True,
        hide_index=True,
        column_config={
            "管子名称": st.column_config.TextColumn(width="medium"),
            "描述": st.column_config.TextColumn(width="large"),
            **{ab: st.column_config.TextColumn(width="small") for ab in antibody_names},
            "固定破膜": st.column_config.TextColumn(width="small"),
            "对照类型": st.column_config.TextColumn(width="small")
        }
    )
    
    # 导出选项
    st.markdown("### 导出矩阵")
    
    col1, col2, col3 = st.columns(3)
    
    with col1:
        csv = df.to_csv(index=False, encoding='utf-8-sig')
        st.download_button(
            label="📥 下载CSV",
            data=csv,
            file_name=f"{st.session_state.current_project}_matrix.csv",
            mime="text/csv",
            use_container_width=True
        )
    
    with col2:
        # 创建简化的表格用于可视化
        simple_df = df.copy()
        # 缩短抗体名称显示
        for ab_name in antibody_names:
            if ab_name in simple_df.columns:
                simple_df[ab_name] = simple_df[ab_name].replace({"✓": "●", "○": ""})
        
        st.dataframe(
            simple_df[["管子名称"] + antibody_names + ["固定破膜", "对照类型"]],
            use_container_width=True,
            hide_index=True
        )
    
    with col3:
        # 创建简单的热图
        st.markdown("### 矩阵热图预览")
        
        # 准备热图数据
        heatmap_data = []
        for tube_name in tube_names:
            tube = st.session_state.tubes[tube_name]
            row = []
            for ab_name in antibody_names:
                row.append(1 if ab_name in tube.antibodies else 0)
            heatmap_data.append(row)
        
        heatmap_df = pd.DataFrame(
            heatmap_data,
            index=tube_names,
            columns=[ab.split('(')[0].strip()[:15] + '...' if len(ab) > 15 else ab.split('(')[0].strip() 
                    for ab in antibody_names]
        )
        
        # 显示为样式化的表格
        def color_cells(val):
            return 'background-color: #10B981' if val == 1 else 'background-color: #F3F4F6'
        
        st.dataframe(
            heatmap_df.style.applymap(color_cells),
            use_container_width=True
        )

def render_mastermix_calculator():
    """母液计算器页面"""
    st.markdown('<div class="section-header">🧪 母液配方计算器</div>', unsafe_allow_html=True)
    
    if not st.session_state.tubes:
        st.warning("请先配置管子")
        return
    
    # 参数设置
    st.markdown("### 实验参数设置")
    
    col1, col2, col3, col4 = st.columns(4)
    
    with col1:
        cell_count = st.number_input(
            "细胞数 (×10⁶)",
            min_value=0.1,
            max_value=100.0,
            value=st.session_state.volumes['cell_count'],
            step=0.5
        )
    
    with col2:
        per_tube = st.number_input(
            "每管体积 (μL)",
            min_value=10.0,
            max_value=200.0,
            value=st.session_state.volumes['per_tube'],
            step=10.0
        )
    
    with col3:
        intracel_volume = st.number_input(
            "胞内染色体积 (μL)",
            min_value=20.0,
            max_value=100.0,
            value=st.session_state.volumes['intracellular_per_tube'],
            step=5.0
        )
    
    with col4:
        extra_tubes = st.number_input(
            "安全余量管数",
            min_value=0,
            max_value=10,
            value=st.session_state.volumes['extra_tubes'],
            step=1
        )
    
    st.session_state.volumes.update({
        'cell_count': cell_count,
        'per_tube': per_tube,
        'intracellular_per_tube': intracel_volume,
        'extra_tubes': extra_tubes
    })
    
    # 统计管子
    surface_tubes = []
    intracellular_tubes = []
    
    for tube_name, tube in st.session_state.tubes.items():
        if tube.needs_fixation:
            intracellular_tubes.append(tube_name)
        elif tube.antibodies:
            surface_tubes.append(tube_name)
    
    # 计算结果
    st.markdown("### 📊 计算结果")
    
    if surface_tubes:
        st.markdown(f"#### 🔬 表面染色母液 (用于 {len(surface_tubes)} 管)")
        
        total_surface_tubes = len(surface_tubes) + extra_tubes
        total_surface_volume = per_tube * total_surface_tubes
        
        surface_data = []
        
        for tube_name in surface_tubes:
            tube = st.session_state.tubes[tube_name]
            for ab_name in tube.antibodies:
                if ab_name in st.session_state.antibodies:
                    ab = st.session_state.antibodies[ab_name]
                    if ab.type in [AntibodyType.SURFACE, AntibodyType.VIABILITY, AntibodyType.FC_BLOCK]:
                        per_tube_vol = (ab.recommended_use * cell_count) / ab.concentration
                        total_vol = per_tube_vol * total_surface_tubes
                        
                        surface_data.append({
                            "抗体": ab.short_name,
                            "荧光": ab.fluorochrome,
                            "每管用量 (μL)": round(per_tube_vol, 2),
                            "总需用量 (μL)": round(total_vol, 2)
                        })
        
        if surface_data:
            surface_df = pd.DataFrame(surface_data).drop_duplicates(subset=["抗体"])
            
            col_surf1, col_surf2 = st.columns([2, 1])
            
            with col_surf1:
                st.dataframe(
                    surface_df,
                    use_container_width=True,
                    hide_index=True
                )
            
            with col_surf2:
                st.info(f"""
                **配制说明:**
                
                - **总体积:** {total_surface_volume} μL
                - **适用管子:** {', '.join(surface_tubes)}
                - **配制步骤:**
                  1. 取洁净EP管
                  2. 按上表加入各抗体
                  3. 用流式染色缓冲液补至{total_surface_volume} μL
                  4. 涡旋混匀，4°C避光保存
                """)
    
    if intracellular_tubes:
        st.markdown(f"#### 🧫 胞内染色工作液 (用于 {len(intracellular_tubes)} 管)")
        
        total_intracel_tubes = len(intracellular_tubes) + extra_tubes
        total_intracel_volume = intracel_volume * total_intracel_tubes
        
        intracel_data = []
        
        for tube_name in intracellular_tubes:
            tube = st.session_state.tubes[tube_name]
            for ab_name in tube.antibodies:
                if ab_name in st.session_state.antibodies:
                    ab = st.session_state.antibodies[ab_name]
                    if ab.type == AntibodyType.INTRACELLULAR:
                        per_tube_vol = (ab.recommended_use * cell_count) / ab.concentration
                        total_vol = per_tube_vol * total_intracel_tubes
                        
                        intracel_data.append({
                            "抗体": ab.short_name,
                            "荧光": ab.fluorochrome,
                            "每管用量 (μL)": round(per_tube_vol, 2),
                            "总需用量 (μL)": round(total_vol, 2)
                        })
        
        if intracel_data:
            intracel_df = pd.DataFrame(intracel_data).drop_duplicates(subset=["抗体"])
            
            col_int1, col_int2 = st.columns([2, 1])
            
            with col_int1:
                st.dataframe(
                    intracel_df,
                    use_container_width=True,
                    hide_index=True
                )
            
            with col_int2:
                st.info(f"""
                **配制说明:**
                
                - **总体积:** {total_intracel_volume} μL
                - **稀释剂:** 1X破膜缓冲液
                - **适用管子:** {', '.join(intracellular_tubes)}
                - **配制步骤:**
                  1. 用1X破膜缓冲液配制
                  2. 按上表加入各抗体
                  3. 用破膜缓冲液补至{total_intracel_volume} μL
                  4. 涡旋混匀，4°C避光保存
                """)

def render_experiment_planner():
    """实验计划页面"""
    st.markdown('<div class="section-header">📋 实验计划生成器</div>', unsafe_allow_html=True)
    
    if not st.session_state.tubes:
        st.warning("请先配置管子")
        return
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.markdown("### 实验组设置")
        
        default_groups = ["Control", "Model", "Treatment"]
        groups_input = st.text_area(
            "实验组名称 (用逗号分隔)",
            value=", ".join(default_groups),
            height=100
        )
        
        groups = [g.strip() for g in groups_input.split(',') if g.strip()]
        
        replicates = st.number_input(
            "每组重复数",
            min_value=1,
            max_value=10,
            value=3,
            step=1
        )
        
        if st.button("🎯 生成实验计划", use_container_width=True):
            st.session_state.experiment_groups = groups
            st.session_state.experiment_replicates = replicates
    
    with col2:
        if 'experiment_groups' in st.session_state:
            st.markdown("### 实验计划概览")
            
            total_samples = len(st.session_state.experiment_groups) * \
                          st.session_state.experiment_replicates * \
                          len(st.session_state.tubes)
            
            st.metric("实验组数", len(st.session_state.experiment_groups))
            st.metric("每组重复", st.session_state.experiment_replicates)
            st.metric("总样品数", total_samples)
    
    if 'experiment_groups' in st.session_state:
        st.markdown("### 📊 详细实验计划")
        
        # 生成计划数据
        plan_data = []
        sample_id = 1
        
        for group in st.session_state.experiment_groups:
            for rep in range(1, st.session_state.experiment_replicates + 1):
                for tube_name, tube in st.session_state.tubes.items():
                    plan_data.append({
                        "样品ID": f"{group[:3]}_R{rep}_{tube_name[:8]}",
                        "实验组": group,
                        "重复": rep,
                        "管子类型": tube_name,
                        "描述": tube.description,
                        "抗体数": len(tube.antibodies),
                        "固定破膜": "是" if tube.needs_fixation else "否",
                        "对照类型": tube.control_type if tube.is_control else "实验管"
                    })
        
        plan_df = pd.DataFrame(plan_data)
        
        # 显示计划表
        st.dataframe(
            plan_df,
            use_container_width=True,
            hide_index=True
        )
        
        # 导出选项
        st.markdown("### 导出选项")
        
        col_exp1, col_exp2 = st.columns(2)
        
        with col_exp1:
            csv = plan_df.to_csv(index=False, encoding='utf-8-sig')
            st.download_button(
                label="📥 下载实验计划 (CSV)",
                data=csv,
                file_name=f"{st.session_state.current_project}_plan.csv",
                mime="text/csv",
                use_container_width=True
            )
        
        with col_exp2:
            # 工作单
            worksheet_df = plan_df[["样品ID", "实验组", "重复", "管子类型", "固定破膜"]].copy()
            worksheet_csv = worksheet_df.to_csv(index=False, encoding='utf-8-sig')
            st.download_button(
                label="📝 下载工作单",
                data=worksheet_csv,
                file_name=f"{st.session_state.current_project}_worksheet.csv",
                mime="text/csv",
                use_container_width=True
            )

def render_protocol():
    """实验方案页面"""
    st.markdown('<div class="section-header">📖 实验方案生成</div>', unsafe_allow_html=True)
    
    if not st.session_state.tubes:
        st.warning("请先配置管子")
        return
    
    # 生成实验方案
    protocol = f"""
# 流式细胞术染色实验方案

## 项目信息
- **项目名称**: {st.session_state.current_project}
- **生成时间**: {datetime.now().strftime('%Y-%m-%d %H:%M')}
- **抗体种类**: {len(st.session_state.antibodies)}种
- **管子配置**: {len(st.session_state.tubes)}种

## 实验步骤

### 1. 样本准备
1. 制备单细胞悬液，用70μm细胞筛过滤
2. 细胞计数，调整浓度至1×10⁷ cells/mL
3. 按实验计划分装细胞到标记好的流式管中
4. 每管分装100μL细胞悬液（约1×10⁶ cells）

### 2. Fc受体阻断与表面染色
1. 配制表面染色母液
2. 向对应管子中加入{st.session_state.volumes['per_tube']}μL母液
3. 4°C避光孵育30分钟
4. 加入1mL预冷染色缓冲液，300g 4°C离心5分钟
5. 弃上清，重复洗涤一次

### 3. 固定与破膜（仅需胞内染色的管子）
1. 每管加入100μL固定液
2. 室温避光孵育20分钟
3. 每管加入1mL 1X破膜缓冲液，300g 4°C离心5分钟
4. 弃上清，重复洗涤一次

### 4. 胞内染色
1. 用1X破膜缓冲液配制胞内抗体工作液
2. 向对应管子中加入{st.session_state.volumes['intracellular_per_tube']}μL工作液
3. 4°C避光孵育45分钟
4. 用1X破膜缓冲液洗涤2次

### 5. 上机检测
1. 所有管子用300μL流式染色缓冲液重悬
2. 过35μm细胞筛网
3. 按上机顺序进行检测

## 质量控制
### 必须设置的对照
1. **未染色对照（Blank）**: 调节电压
2. **单阳对照**: 用于荧光补偿
3. **FMO对照**: 用于准确设门

### 注意事项
- 所有操作避光进行
- 离心条件: 300g, 4°C, 5分钟
- 抗体现配现用
- 设置正确的补偿
"""
    
    # 显示方案
    st.markdown(protocol)
    
    # 导出方案
    st.markdown("### 导出实验方案")
    
    st.download_button(
        label="📄 下载实验方案 (TXT)",
        data=protocol,
        file_name=f"{st.session_state.current_project}_protocol.txt",
        mime="text/plain",
        use_container_width=True
    )

def render_sidebar():
    """侧边栏"""
    with st.sidebar:
        st.markdown("## 🔬 流式染色矩阵设计器")
        st.markdown("---")
        
        # 项目信息
        st.markdown(f"**项目:** {st.session_state.current_project}")
        st.markdown(f"**抗体:** {len(st.session_state.antibodies)}种")
        st.markdown(f"**管子:** {len(st.session_state.tubes)}种")
        
        st.markdown("---")
        
        # 导航菜单
        st.markdown("### 📋 导航菜单")
        
        pages = {
            "🏠 仪表盘": "dashboard",
            "🧪 抗体管理": "antibodies",
            "🧫 管子设计": "tubes",
            "🔢 染色矩阵": "matrix",
            "🧪 母液计算": "mastermix",
            "📋 实验计划": "planner",
            "📖 实验方案": "protocol"
        }
        
        if 'page' not in st.session_state:
            st.session_state.page = "dashboard"
        
        selected = st.radio(
            "选择页面",
            options=list(pages.keys()),
            index=list(pages.values()).index(st.session_state.page) if st.session_state.page in pages.values() else 0,
            label_visibility="collapsed"
        )
        
        st.session_state.page = pages[selected]
        
        st.markdown("---")
        
        # 快速操作
        st.markdown("### ⚡ 快速操作")
        
        if st.button("🔄 重置项目", use_container_width=True, type="secondary"):
            for key in ['antibodies', 'tubes']:
                if key in st.session_state:
                    st.session_state[key] = {}
            st.session_state.page = "dashboard"
            st.rerun()
        
        st.markdown("---")
        st.markdown("### ℹ️ 使用说明")
        st.markdown("""
        1. 从"抗体管理"开始添加抗体
        2. 在"管子设计"中配置实验管
        3. 查看"染色矩阵"确认配置
        4. 使用"母液计算"获取配方
        5. 生成"实验计划"和"实验方案"
        """)

def main():
    """主函数"""
    # 初始化session state
    init_session()
    
    # 渲染侧边栏
    render_sidebar()
    
    # 根据选择渲染页面
    page = st.session_state.get('page', 'dashboard')
    
    if page == "dashboard":
        render_dashboard()
    elif page == "antibodies":
        render_antibody_management()
    elif page == "tubes":
        render_tube_design()
    elif page == "matrix":
        render_matrix()
    elif page == "mastermix":
        render_mastermix_calculator()
    elif page == "planner":
        render_experiment_planner()
    elif page == "protocol":
        render_protocol()

if __name__ == "__main__":
    main()
