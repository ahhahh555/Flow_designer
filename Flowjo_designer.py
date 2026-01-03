# app.py
import streamlit as st
import pandas as pd
import numpy as np
from typing import Dict, List, Optional, Tuple
import json
from dataclasses import dataclass, field, asdict
from enum import Enum
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import base64
from datetime import datetime
import io

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
    .matrix-cell {
        text-align: center;
        font-weight: bold;
    }
    .matrix-cell-present {
        background-color: #10B981 !important;
        color: white !important;
    }
    .matrix-cell-absent {
        background-color: #F3F4F6 !important;
        color: #6B7280 !important;
    }
    .antibody-card {
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
        color: white;
        padding: 1rem;
        border-radius: 12px;
        margin: 0.5rem 0;
    }
    .tube-card {
        background: linear-gradient(135deg, #f093fb 0%, #f5576c 100%);
        color: white;
        padding: 1rem;
        border-radius: 12px;
        margin: 0.5rem 0;
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
    color: str = "#3B82F6"  # 可视化颜色
    
    def __post_init__(self):
        if not self.short_name:
            # 自动生成简称
            if "Anti-" in self.name:
                self.short_name = self.name.split("Anti-")[-1].split()[0]
            elif " " in self.name:
                self.short_name = ''.join([word[0] for word in self.name.split() if word[0].isupper()])
            else:
                self.short_name = self.name[:10]

@dataclass
class TubeConfiguration:
    """管子配置"""
    name: str
    description: str
    antibodies: List[str] = field(default_factory=list)
    needs_fixation: bool = False
    is_control: bool = False
    control_type: str = ""  # FMO, Isotype, Single, Blank
    color: str = "#10B981"  # 可视化颜色

class FlowMatrixDesigner:
    """流式矩阵设计器"""
    
    def __init__(self):
        if 'antibodies' not in st.session_state:
            st.session_state.antibodies = {}
        if 'tubes' not in st.session_state:
            st.session_state.tubes = {}
        if 'current_project' not in st.session_state:
            st.session_state.current_project = f"Flow_Project_{datetime.now().strftime('%Y%m%d_%H%M')}"
        if 'volumes' not in st.session_state:
            st.session_state.volumes = {
                'per_tube': 100.0,
                'intracellular_per_tube': 50.0,
                'cell_count': 1.0,  # 10⁶ cells
                'extra_tubes': 2
            }
    
    def load_standard_antibodies(self):
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
                notes="Fc受体阻断剂，所有管子都需要添加",
                color="#FF6B6B"
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
                notes="死活染料，建议1:1000稀释使用",
                color="#4ECDC4"
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
                notes="白细胞标记，BB515通道（类似FITC）",
                color="#45B7D1"
            ),
            Antibody(
                name="α-SMA AF647",
                short_name="α-SMA",
                fluorochrome="AF647",
                target="α-Smooth Muscle Actin",
                clone="1A4",
                concentration=200.0,
                recommended_use=0.5,
                type=AntibodyType.INTRACELLULAR,
                catalog_number="561847",
                lot_number="345678",
                notes="胞内染色，需固定破膜，AF647通道（类似APC）",
                color="#96CEB4"
            )
        ]
        
        for ab in standard_antibodies:
            st.session_state.antibodies[ab.name] = ab
        
        return True
    
    def load_standard_tubes(self):
        """加载标准管子配置"""
        standard_tubes = {
            "Blank": TubeConfiguration(
                name="Blank",
                description="未染色对照，用于调节电压和检测自发荧光",
                is_control=True,
                control_type="Blank",
                color="#FFD166"
            ),
            "FcX_Only": TubeConfiguration(
                name="FcX_Only",
                description="仅Fc阻断对照",
                antibodies=["TruStain FcX™ (anti-mouse CD16/32) Antibody"],
                is_control=True,
                control_type="Single",
                color="#EF476F"
            ),
            "Live_Only": TubeConfiguration(
                name="Live_Only",
                description="死活染料单阳（用于补偿调节）",
                antibodies=["TruStain FcX™ (anti-mouse CD16/32) Antibody", "Live/Dye eF780"],
                is_control=True,
                control_type="Single",
                color="#118AB2"
            ),
            "CD45_Only": TubeConfiguration(
                name="CD45_Only",
                description="CD45单阳（用于补偿调节）",
                antibodies=["TruStain FcX™ (anti-mouse CD16/32) Antibody", "BB515 Rat Anti-Mouse CD45"],
                is_control=True,
                control_type="Single",
                color="#06D6A0"
            ),
            "αSMA_Only": TubeConfiguration(
                name="αSMA_Only",
                description="α-SMA单阳（需破膜，用于补偿调节）",
                antibodies=["TruStain FcX™ (anti-mouse CD16/32) Antibody", "α-SMA AF647"],
                needs_fixation=True,
                is_control=True,
                control_type="Single",
                color="#073B4C"
            ),
            "FMO_αSMA": TubeConfiguration(
                name="FMO_αSMA",
                description="荧光减一对照（不含α-SMA，用于准确设门）",
                antibodies=["TruStain FcX™ (anti-mouse CD16/32) Antibody", "Live/Dye eF780", "BB515 Rat Anti-Mouse CD45"],
                is_control=True,
                control_type="FMO",
                color="#7209B7"
            ),
            "Full_Stain": TubeConfiguration(
                name="Full_Stain",
                description="全染实验管（包含所有抗体）",
                antibodies=["TruStain FcX™ (anti-mouse CD16/32) Antibody", "Live/Dye eF780", "BB515 Rat Anti-Mouse CD45", "α-SMA AF647"],
                needs_fixation=True,
                is_control=False,
                color="#F72585"
            )
        }
        
        st.session_state.tubes = standard_tubes
        return True

def init_session():
    """初始化session state"""
    designer = FlowMatrixDesigner()
    
    # 如果没有数据，加载标准配置
    if not st.session_state.antibodies:
        if st.sidebar.button("📥 加载标准实验模板", use_container_width=True):
            designer.load_standard_antibodies()
            designer.load_standard_tubes()
            st.rerun()

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
        st.markdown(f"""
        <div style="
            background: linear-gradient(135deg, {color}20 0%, {color}40 100%);
            border-left: 5px solid {color};
            padding: 1rem;
            border-radius: 8px;
            margin: 0.5rem 0;
        ">
            <div style="display: flex; justify-content: space-between; align-items: center;">
                <div>
                    <h4 style="margin: 0; color: {color};">{antibody.name}</h4>
                    <p style="margin: 0.2rem 0; color: #6B7280; font-size: 0.9rem;">
                        <strong>靶标:</strong> {antibody.target} | 
                        <strong>荧光:</strong> {antibody.fluorochrome} | 
                        <strong>克隆:</strong> {antibody.clone}
                    </p>
                </div>
                <span style="
                    background-color: {color};
                    color: white;
                    padding: 0.2rem 0.8rem;
                    border-radius: 12px;
                    font-size: 0.8rem;
                ">{antibody.type.value}</span>
            </div>
            <div style="margin-top: 0.8rem; font-size: 0.9rem;">
                <p style="margin: 0.3rem 0;">
                    <strong>浓度:</strong> {antibody.concentration} μg/mL | 
                    <strong>用量:</strong> {antibody.recommended_use} μg/10⁶ cells
                </p>
                <p style="margin: 0.3rem 0;">
                    <strong>货号:</strong> {antibody.catalog_number} | 
                    <strong>批号:</strong> {antibody.lot_number}
                </p>
                {f'<p style="margin: 0.3rem 0;"><strong>备注:</strong> {antibody.notes}</p>' if antibody.notes else ''}
            </div>
        </div>
        """, unsafe_allow_html=True)

def display_tube_card(tube: TubeConfiguration):
    """显示管子卡片"""
    control_colors = {
        "FMO": "#8B5CF6",
        "Isotype": "#F59E0B",
        "Single": "#10B981",
        "Blank": "#6B7280",
        "": "#3B82F6"
    }
    
    bg_color = control_colors.get(tube.control_type, "#3B82F6") if tube.is_control else "#3B82F6"
    
    with st.container():
        st.markdown(f"""
        <div style="
            background: linear-gradient(135deg, {bg_color}20 0%, {bg_color}40 100%);
            border-left: 5px solid {bg_color};
            padding: 1rem;
            border-radius: 8px;
            margin: 0.5rem 0;
        ">
            <div style="display: flex; justify-content: space-between; align-items: center;">
                <div>
                    <h4 style="margin: 0; color: {bg_color};">{tube.name}</h4>
                    <p style="margin: 0.2rem 0; color: #6B7280; font-size: 0.9rem;">
                        {tube.description}
                    </p>
                </div>
                <div>
                    {f'<span style="background-color: {bg_color}; color: white; padding: 0.2rem 0.8rem; border-radius: 12px; font-size: 0.8rem; margin-left: 0.5rem;">{tube.control_type}对照</span>' if tube.is_control else ''}
                    {f'<span style="background-color: #F59E0B; color: white; padding: 0.2rem 0.8rem; border-radius: 12px; font-size: 0.8rem; margin-left: 0.5rem;">需破膜</span>' if tube.needs_fixation else ''}
                </div>
            </div>
            <div style="margin-top: 0.8rem;">
                <p style="margin: 0.3rem 0; font-size: 0.9rem;">
                    <strong>抗体 ({len(tube.antibodies)}):</strong> {', '.join([ab.split('(')[0].strip() for ab in tube.antibodies]) if tube.antibodies else '无'}
                </p>
            </div>
        </div>
        """, unsafe_allow_html=True)

def render_antibody_management():
    """抗体管理页面"""
    st.markdown('<div class="section-header">🧪 抗体库管理</div>', unsafe_allow_html=True)
    
    col1, col2 = st.columns([1, 2])
    
    with col1:
        st.markdown("### 添加新抗体")
        
        with st.form("add_antibody_form"):
            name = st.text_input("抗体全名*", placeholder="例如: TruStain FcX™ (anti-mouse CD16/32) Antibody")
            short_name = st.text_input("简称", placeholder="例如: FcX")
            fluorochrome = st.text_input("荧光染料", placeholder="例如: AF647, BB515, FITC")
            target = st.text_input("靶标", placeholder="例如: CD45, α-SMA")
            clone = st.text_input("克隆号", placeholder="例如: 30-F11, 1A4")
            
            col_a, col_b = st.columns(2)
            with col_a:
                concentration = st.number_input("浓度 (μg/mL)*", min_value=0.0, value=200.0, step=10.0)
            with col_b:
                recommended_use = st.number_input("用量 (μg/10⁶ cells)*", min_value=0.0, value=0.5, step=0.1)
            
            antibody_type = st.selectbox(
                "抗体类型*",
                options=[t.value for t in AntibodyType],
                index=0
            )
            
            catalog_number = st.text_input("货号", placeholder="例如: 101320")
            lot_number = st.text_input("批号", placeholder="例如: B123456")
            notes = st.text_area("备注", placeholder="例如: 建议1:100稀释使用")
            
            submitted = st.form_submit_button("✅ 添加抗体", use_container_width=True)
            
            if submitted:
                if not name:
                    st.error("抗体名称不能为空！")
                else:
                    # 转换类型字符串为枚举
                    type_map = {t.value: t for t in AntibodyType}
                    
                    antibody = Antibody(
                        name=name,
                        short_name=short_name,
                        fluorochrome=fluorochrome,
                        target=target,
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
                    st.rerun()
    
    with col2:
        st.markdown(f"### 抗体库 ({len(st.session_state.antibodies)}种)")
        
        if not st.session_state.antibodies:
            st.info("暂无抗体数据，请先添加抗体或加载标准模板")
        else:
            search_term = st.text_input("🔍 搜索抗体", placeholder="输入名称、靶标或荧光染料搜索")
            
            filtered_antibodies = {}
            if search_term:
                for name, ab in st.session_state.antibodies.items():
                    if (search_term.lower() in name.lower() or 
                        search_term.lower() in ab.target.lower() or 
                        search_term.lower() in ab.fluorochrome.lower()):
                        filtered_antibodies[name] = ab
            else:
                filtered_antibodies = st.session_state.antibodies
            
            if filtered_antibodies:
                for ab in filtered_antibodies.values():
                    display_antibody_card(ab)
                    
                    col_del, _ = st.columns([1, 5])
                    with col_del:
                        if st.button(f"删除 {ab.short_name}", key=f"del_{ab.name}", type="secondary", use_container_width=True):
                            del st.session_state.antibodies[ab.name]
                            st.rerun()
            else:
                st.warning("未找到匹配的抗体")

def render_tube_design():
    """管子设计页面"""
    st.markdown('<div class="section-header">🧫 实验管子设计</div>', unsafe_allow_html=True)
    
    col1, col2 = st.columns([1, 2])
    
    with col1:
        st.markdown("### 创建新管子配置")
        
        with st.form("add_tube_form"):
            tube_name = st.text_input("管子名称*", placeholder="例如: Full_Stain, FMO_αSMA")
            description = st.text_area("描述", placeholder="例如: 全染实验管，包含所有抗体")
            
            col_type1, col_type2 = st.columns(2)
            with col_type1:
                is_control = st.checkbox("是否为对照管")
            with col_type2:
                needs_fixation = st.checkbox("需要固定破膜")
            
            control_type = ""
            if is_control:
                control_type = st.selectbox(
                    "对照类型",
                    options=["FMO", "Isotype", "Single", "Blank"],
                    index=0
                )
            
            # 选择抗体
            if st.session_state.antibodies:
                st.markdown("### 选择抗体")
                selected_antibodies = []
                
                for name, ab in st.session_state.antibodies.items():
                    if st.checkbox(f"{ab.name} ({ab.fluorochrome})", key=f"tube_{tube_name}_{name}"):
                        selected_antibodies.append(name)
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
                    st.rerun()
    
    with col2:
        st.markdown(f"### 管子配置 ({len(st.session_state.tubes)}种)")
        
        if not st.session_state.tubes:
            st.info("暂无管子配置，请先创建管子或加载标准模板")
        else:
            for tube in st.session_state.tubes.values():
                display_tube_card(tube)
                
                col_edit, col_del = st.columns(2)
                with col_edit:
                    if st.button(f"编辑 {tube.name}", key=f"edit_{tube.name}", type="secondary", use_container_width=True):
                        st.session_state.editing_tube = tube.name
                        st.rerun()
                with col_del:
                    if st.button(f"删除 {tube.name}", key=f"del_tube_{tube.name}", type="secondary", use_container_width=True):
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
            row[ab_name] = "✓" if ab_name in tube.antibodies else ""
        
        row["固定破膜"] = "是" if tube.needs_fixation else "否"
        row["对照类型"] = tube.control_type if tube.is_control else "实验管"
        matrix_data.append(row)
    
    # 创建DataFrame
    df = pd.DataFrame(matrix_data)
    
    # 显示矩阵
    col1, col2 = st.columns([3, 1])
    
    with col1:
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
    
    with col2:
        st.markdown("### 图例说明")
        st.markdown("""
        <div style="background-color: #f0f8ff; padding: 1rem; border-radius: 8px;">
            <p><span style="color: #10B981; font-weight: bold;">✓</span> 表示该抗体存在于管子中</p>
            <p><span style="color: #6B7280; font-weight: bold;"></span> 表示该抗体不存在</p>
            <p><span style="color: #F59E0B; font-weight: bold;">需破膜</span> 需要固定破膜步骤</p>
            <p><strong>对照类型:</strong></p>
            <ul>
                <li>FMO: 荧光减一对照</li>
                <li>Single: 单阳对照</li>
                <li>Isotype: 同型对照</li>
                <li>Blank: 空白对照</li>
            </ul>
        </div>
        """, unsafe_allow_html=True)
    
    # 可视化矩阵
    st.markdown("### 矩阵可视化")
    
    # 准备热图数据
    heatmap_data = []
    for tube_name in tube_names:
        tube = st.session_state.tubes[tube_name]
        row = []
        for ab_name in antibody_names:
            row.append(1 if ab_name in tube.antibodies else 0)
        heatmap_data.append(row)
    
    heatmap_data = np.array(heatmap_data)
    
    # 创建热图
    fig = go.Figure(data=go.Heatmap(
        z=heatmap_data,
        x=[ab.split('(')[0].strip()[:20] + '...' if len(ab) > 20 else ab.split('(')[0].strip() for ab in antibody_names],
        y=tube_names,
        colorscale=[[0, '#F3F4F6'], [1, '#10B981']],
        showscale=False,
        text=[[f"{tube_names[i]}\n{antibody_names[j]}\n{'存在' if heatmap_data[i][j] else '不存在'}" 
               for j in range(len(antibody_names))] for i in range(len(tube_names))],
        hoverinfo="text"
    ))
    
    fig.update_layout(
        title="染色矩阵热图",
        xaxis_title="抗体",
        yaxis_title="管子类型",
        height=400,
        margin=dict(l=20, r=20, t=40, b=20)
    )
    
    st.plotly_chart(fig, use_container_width=True)
    
    # 导出选项
    st.markdown("### 导出矩阵")
    col_exp1, col_exp2, col_exp3 = st.columns(3)
    
    with col_exp1:
        if st.button("📥 导出为CSV", use_container_width=True):
            csv = df.to_csv(index=False, encoding='utf-8-sig')
            st.download_button(
                label="下载CSV文件",
                data=csv,
                file_name=f"{st.session_state.current_project}_matrix.csv",
                mime="text/csv"
            )
    
    with col_exp2:
        if st.button("📊 导出为Excel", use_container_width=True):
            buffer = io.BytesIO()
            with pd.ExcelWriter(buffer, engine='openpyxl') as writer:
                df.to_excel(writer, index=False, sheet_name='染色矩阵')
            st.download_button(
                label="下载Excel文件",
                data=buffer.getvalue(),
                file_name=f"{st.session_state.current_project}_matrix.xlsx",
                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
            )
    
    with col_exp3:
        if st.button("📋 复制到剪贴板", use_container_width=True):
            st.code(df.to_string(index=False), language="text")

def render_mastermix_calculator():
    """母液计算器页面"""
    st.markdown('<div class="section-header">🧪 母液配方计算器</div>', unsafe_allow_html=True)
    
    if not st.session_state.tubes:
        st.warning("请先配置管子")
        return
    
    # 参数设置
    st.markdown("### 实验参数设置")
    
    col_param1, col_param2, col_param3, col_param4 = st.columns(4)
    
    with col_param1:
        cell_count = st.number_input(
            "细胞数 (×10⁶)",
            min_value=0.1,
            max_value=100.0,
            value=st.session_state.volumes['cell_count'],
            step=0.5,
            help="每管的细胞数量"
        )
    
    with col_param2:
        per_tube = st.number_input(
            "每管体积 (μL)",
            min_value=10.0,
            max_value=200.0,
            value=st.session_state.volumes['per_tube'],
            step=10.0,
            help="每管加入的抗体混合液体积"
        )
    
    with col_param3:
        intracel_volume = st.number_input(
            "胞内染色体积 (μL)",
            min_value=20.0,
            max_value=100.0,
            value=st.session_state.volumes['intracellular_per_tube'],
            step=5.0,
            help="每管胞内染色抗体工作液体积"
        )
    
    with col_param4:
        extra_tubes = st.number_input(
            "安全余量管数",
            min_value=0,
            max_value=10,
            value=st.session_state.volumes['extra_tubes'],
            step=1,
            help="多配的管数，用于弥补损耗"
        )
    
    st.session_state.volumes.update({
        'cell_count': cell_count,
        'per_tube': per_tube,
        'intracellular_per_tube': intracel_volume,
        'extra_tubes': extra_tubes
    })
    
    # 统计需要不同混合液的管子
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
        
        surface_results = []
        
        for tube_name in surface_tubes:
            tube = st.session_state.tubes[tube_name]
            for ab_name in tube.antibodies:
                if ab_name in st.session_state.antibodies:
                    ab = st.session_state.antibodies[ab_name]
                    if ab.type in [AntibodyType.SURFACE, AntibodyType.VIABILITY, AntibodyType.FC_BLOCK]:
                        per_tube_vol = (ab.recommended_use * cell_count) / ab.concentration
                        total_vol = per_tube_vol * total_surface_tubes
                        
                        surface_results.append({
                            "抗体": ab.short_name,
                            "浓度": f"{ab.concentration} μg/mL",
                            "每管用量": f"{per_tube_vol:.2f} μL",
                            "总需用量": f"{total_vol:.2f} μL",
                            "备注": ab.notes if hasattr(ab, 'notes') else ""
                        })
        
        if surface_results:
            surface_df = pd.DataFrame(surface_results).drop_duplicates(subset=["抗体"])
            
            col_surf1, col_surf2 = st.columns([2, 1])
            
            with col_surf1:
                st.dataframe(
                    surface_df,
                    use_container_width=True,
                    hide_index=True
                )
            
            with col_surf2:
                st.markdown(f"""
                <div class="info-box">
                    <h4>配制说明</h4>
                    <p><strong>总体积:</strong> {total_surface_volume} μL</p>
                    <p><strong>配制步骤:</strong></p>
                    <ol>
                        <li>取洁净EP管</li>
                        <li>按上表加入各抗体</li>
                        <li>用流式染色缓冲液补至{total_surface_volume} μL</li>
                        <li>涡旋混匀，4°C避光保存</li>
                    </ol>
                    <p><strong>适用管子:</strong> {', '.join(surface_tubes)}</p>
                </div>
                """, unsafe_allow_html=True)
    
    if intracellular_tubes:
        st.markdown(f"#### 🧫 胞内染色工作液 (用于 {len(intracellular_tubes)} 管)")
        
        total_intracel_tubes = len(intracellular_tubes) + extra_tubes
        total_intracel_volume = intracel_volume * total_intracel_tubes
        
        intracel_results = []
        
        for tube_name in intracellular_tubes:
            tube = st.session_state.tubes[tube_name]
            for ab_name in tube.antibodies:
                if ab_name in st.session_state.antibodies:
                    ab = st.session_state.antibodies[ab_name]
                    if ab.type == AntibodyType.INTRACELLULAR:
                        per_tube_vol = (ab.recommended_use * cell_count) / ab.concentration
                        total_vol = per_tube_vol * total_intracel_tubes
                        
                        intracel_results.append({
                            "抗体": ab.short_name,
                            "浓度": f"{ab.concentration} μg/mL",
                            "每管用量": f"{per_tube_vol:.2f} μL",
                            "总需用量": f"{total_vol:.2f} μL",
                            "备注": ab.notes if hasattr(ab, 'notes') else ""
                        })
        
        if intracel_results:
            intracel_df = pd.DataFrame(intracel_results).drop_duplicates(subset=["抗体"])
            
            col_int1, col_int2 = st.columns([2, 1])
            
            with col_int1:
                st.dataframe(
                    intracel_df,
                    use_container_width=True,
                    hide_index=True
                )
            
            with col_int2:
                st.markdown(f"""
                <div class="info-box">
                    <h4>配制说明</h4>
                    <p><strong>总体积:</strong> {total_intracel_volume} μL</p>
                    <p><strong>稀释剂:</strong> 1X破膜缓冲液</p>
                    <p><strong>配制步骤:</strong></p>
                    <ol>
                        <li>用1X破膜缓冲液配制</li>
                        <li>按上表加入各抗体</li>
                        <li>用破膜缓冲液补至{total_intracel_volume} μL</li>
                        <li>涡旋混匀，4°C避光保存</li>
                        <li>建议使用前过滤</li>
                    </ol>
                    <p><strong>适用管子:</strong> {', '.join(intracellular_tubes)}</p>
                </div>
                """, unsafe_allow_html=True)
    
    # 抗体用量统计图表
    if st.session_state.antibodies:
        st.markdown("### 📈 抗体用量统计")
        
        usage_data = []
        for ab_name, ab in st.session_state.antibodies.items():
            # 计算总使用管数
            usage_count = 0
            for tube in st.session_state.tubes.values():
                if ab_name in tube.antibodies:
                    if tube.needs_fixation and ab.type == AntibodyType.INTRACELLULAR:
                        usage_count += len(intracellular_tubes) + extra_tubes
                    elif not tube.needs_fixation and ab.type != AntibodyType.INTRACELLULAR:
                        usage_count += len(surface_tubes) + extra_tubes
            
            if usage_count > 0:
                per_tube_vol = (ab.recommended_use * cell_count) / ab.concentration
                total_vol = per_tube_vol * usage_count
                
                usage_data.append({
                    "抗体": ab.short_name,
                    "荧光": ab.fluorochrome,
                    "类型": ab.type.value,
                    "使用管数": usage_count,
                    "总用量 (μL)": total_vol
                })
        
        if usage_data:
            usage_df = pd.DataFrame(usage_data)
            
            # 创建图表
            fig = make_subplots(
                rows=1, cols=2,
                subplot_titles=('抗体使用管数', '抗体总用量 (μL)'),
                specs=[[{"type": "bar"}, {"type": "bar"}]]
            )
            
            # 使用管数条形图
            fig.add_trace(
                go.Bar(
                    x=usage_df["抗体"],
                    y=usage_df["使用管数"],
                    name="使用管数",
                    marker_color='#3B82F6',
                    text=usage_df["使用管数"],
                    textposition='auto'
                ),
                row=1, col=1
            )
            
            # 总用量条形图
            fig.add_trace(
                go.Bar(
                    x=usage_df["抗体"],
                    y=usage_df["总用量 (μL)"],
                    name="总用量",
                    marker_color='#10B981',
                    text=usage_df["总用量 (μL)"].round(2),
                    textposition='auto'
                ),
                row=1, col=2
            )
            
            fig.update_layout(
                height=400,
                showlegend=False,
                margin=dict(l=20, r=20, t=40, b=20)
            )
            
            fig.update_xaxes(tickangle=45)
            
            st.plotly_chart(fig, use_container_width=True)

def render_experiment_planner():
    """实验计划页面"""
    st.markdown('<div class="section-header">📋 实验计划生成器</div>', unsafe_allow_html=True)
    
    if not st.session_state.tubes:
        st.warning("请先配置管子")
        return
    
    col_plan1, col_plan2 = st.columns(2)
    
    with col_plan1:
        st.markdown("### 实验组设置")
        
        default_groups = ["Control", "Model", "Treatment_1", "Treatment_2"]
        groups_input = st.text_area(
            "实验组名称",
            value=", ".join(default_groups),
            help="用逗号分隔多个实验组",
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
        
        randomize = st.checkbox("随机化上机顺序", value=True)
        
        if st.button("🎯 生成实验计划", use_container_width=True):
            st.session_state.groups = groups
            st.session_state.replicates = replicates
    
    with col_plan2:
        if 'groups' in st.session_state:
            st.markdown("### 实验计划概览")
            
            total_samples = len(st.session_state.groups) * replicates * len(st.session_state.tubes)
            
            st.metric("实验组数", len(st.session_state.groups))
            st.metric("每组重复", replicates)
            st.metric("管子类型", len(st.session_state.tubes))
            st.metric("总样品数", total_samples)
    
    if 'groups' in st.session_state:
        st.markdown("### 📊 详细实验计划")
        
        # 生成计划数据
        plan_data = []
        sample_id = 1
        
        for group in st.session_state.groups:
            for rep in range(1, replicates + 1):
                for tube_name, tube in st.session_state.tubes.items():
                    plan_data.append({
                        "样品ID": f"{group[:3]}_R{rep}_{tube_name[:10]}",
                        "实验组": group,
                        "重复": rep,
                        "管子类型": tube_name,
                        "描述": tube.description,
                        "抗体": ", ".join([ab.split('(')[0].strip() for ab in tube.antibodies]),
                        "固定破膜": "是" if tube.needs_fixation else "否",
                        "对照类型": tube.control_type if tube.is_control else "实验管",
                        "上机顺序": sample_id
                    })
                    sample_id += 1
        
        plan_df = pd.DataFrame(plan_data)
        
        # 随机化上机顺序
        if randomize:
            np.random.seed(42)  # 固定随机种子以便重现
            plan_df = plan_df.sample(frac=1).reset_index(drop=True)
            plan_df["上机顺序"] = range(1, len(plan_df) + 1)
        
        # 显示计划表
        st.dataframe(
            plan_df,
            use_container_width=True,
            hide_index=True,
            column_config={
                "样品ID": st.column_config.TextColumn(width="small"),
                "实验组": st.column_config.TextColumn(width="medium"),
                "重复": st.column_config.NumberColumn(width="small"),
                "管子类型": st.column_config.TextColumn(width="medium"),
                "描述": st.column_config.TextColumn(width="large"),
                "抗体": st.column_config.TextColumn(width="large"),
                "固定破膜": st.column_config.TextColumn(width="small"),
                "对照类型": st.column_config.TextColumn(width="small"),
                "上机顺序": st.column_config.NumberColumn(width="small")
            }
        )
        
        # 导出选项
        st.markdown("### 导出选项")
        col_exp1, col_exp2, col_exp3 = st.columns(3)
        
        with col_exp1:
            csv = plan_df.to_csv(index=False, encoding='utf-8-sig')
            st.download_button(
                label="📥 下载实验计划 (CSV)",
                data=csv,
                file_name=f"{st.session_state.current_project}_experiment_plan.csv",
                mime="text/csv",
                use_container_width=True
            )
        
        with col_exp2:
            # 上机顺序表
            run_order_df = plan_df[["样品ID", "实验组", "管子类型", "上机顺序"]].copy()
            run_order_csv = run_order_df.to_csv(index=False, encoding='utf-8-sig')
            st.download_button(
                label="📋 下载上机顺序表",
                data=run_order_csv,
                file_name=f"{st.session_state.current_project}_run_order.csv",
                mime="text/csv",
                use_container_width=True
            )
        
        with col_exp3:
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
        
        # 可视化实验布局
        st.markdown("### 🎨 实验布局可视化")
        
        fig = go.Figure()
        
        colors = px.colors.qualitative.Set3
        
        for i, group in enumerate(st.session_state.groups):
            group_data = plan_df[plan_df["实验组"] == group]
            
            fig.add_trace(go.Bar(
                x=[group],
                y=[len(group_data)],
                name=group,
                marker_color=colors[i % len(colors)],
                text=f"{len(group_data)}个样品",
                textposition='auto'
            ))
        
        fig.update_layout(
            title="各实验组样品数量",
            xaxis_title="实验组",
            yaxis_title="样品数量",
            barmode='group',
            height=400,
            margin=dict(l=20, r=20, t=40, b=20)
        )
        
        st.plotly_chart(fig, use_container_width=True)

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

## 实验设计概览
### 抗体列表
{chr(10).join([f"- {ab.name} ({ab.fluorochrome}): {ab.target}, {ab.concentration} μg/mL" for ab in st.session_state.antibodies.values()])}

### 管子配置
{chr(10).join([f"- {tube.name}: {tube.description}" for tube in st.session_state.tubes.values()])}

## 实验步骤

### 1. 样本准备
1. 制备单细胞悬液，用70μm细胞筛过滤
2. 细胞计数，用预冷流式染色缓冲液调整浓度至1×10⁷ cells/mL
3. 按实验计划分装细胞到标记好的流式管中
4. 每管分装100μL细胞悬液（约1×10⁶ cells）

### 2. Fc受体阻断与表面染色
1. 按"母液配方计算器"的结果配制表面染色母液
2. 向对应管子中加入{st.session_state.volumes['per_tube']}μL表面染色母液
3. 4°C避光孵育30分钟
4. 每管加入1mL预冷染色缓冲液，300g 4°C离心5分钟
5. 弃上清，重复洗涤一次

### 3. 固定与破膜（仅需胞内染色的管子）
1. 每管加入100μL固定液（如BD Cytofix）
2. 室温避光孵育20分钟
3. 每管加入1mL 1X破膜缓冲液，300g 4°C离心5分钟
4. 弃上清，重复洗涤一次

### 4. 胞内染色
1. 用1X破膜缓冲液配制胞内抗体工作液
2. 向对应管子中加入{st.session_state.volumes['intracellular_per_tube']}μL抗体工作液
3. 4°C避光孵育45分钟
4. 用1X破膜缓冲液洗涤2次

### 5. 上机检测
1. 所有管子用300μL流式染色缓冲液重悬
2. 过35μm细胞筛网
3. 按上机顺序表进行检测

## 质量控制
### 必须设置的对照
1. **未染色对照（Blank）**: 调节电压，确定自发荧光水平
2. **单阳对照**: 用于荧光补偿调节
3. **FMO对照**: 用于准确设门
4. **同型对照（可选）**: 检测非特异性结合

### 上机注意事项
1. 上机前充分混匀样本
2. 使用低速采集模式（<500 events/sec）
3. 先调节电压，再设置补偿
4. 每管至少采集10,000个目标细胞事件

## 数据分析建议
### 设门策略
1. FSC-A vs SSC-A: 圈出有核细胞，排除碎片
2. FSC-H vs FSC-A: 排除粘连体
3. 死活染料阴性: 圈出活细胞
4. 根据实验目标分析相应群体

### 数据记录
1. 记录各抗体货号、批号、用量
2. 保存原始FCS文件
3. 记录仪器设置和补偿矩阵

## 备注
- 所有操作避光进行
- 抗体工作液现配现用
- 离心条件: 300g, 4°C, 5分钟
- 实验完成后及时分析数据
"""
    
    # 显示方案
    st.markdown(protocol)
    
    # 导出方案
    st.markdown("### 导出实验方案")
    
    col_prot1, col_prot2 = st.columns(2)
    
    with col_prot1:
        # 导出为文本文件
        st.download_button(
            label="📄 下载实验方案 (TXT)",
            data=protocol,
            file_name=f"{st.session_state.current_project}_protocol.txt",
            mime="text/plain",
            use_container_width=True
        )
    
    with col_prot2:
        # 导出为PDF（通过HTML转换）
        html_protocol = f"""
        <html>
        <head>
            <meta charset="UTF-8">
            <title>流式染色实验方案 - {st.session_state.current_project}</title>
            <style>
                body {{ font-family: Arial, sans-serif; margin: 40px; }}
                h1 {{ color: #1E3A8A; border-bottom: 2px solid #3B82F6; }}
                h2 {{ color: #2563EB; margin-top: 30px; }}
                .section {{ margin-bottom: 20px; }}
                table {{ border-collapse: collapse; width: 100%; margin: 20px 0; }}
                th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
                th {{ background-color: #f2f2f2; }}
            </style>
        </head>
        <body>
            <h1>流式细胞术染色实验方案</h1>
            <p><strong>项目名称:</strong> {st.session_state.current_project}</p>
            <p><strong>生成时间:</strong> {datetime.now().strftime('%Y-%m-%d %H:%M')}</p>
            
            <h2>抗体列表</h2>
            <table>
                <tr><th>抗体名称</th><th>荧光</th><th>靶标</th><th>浓度</th></tr>
                {"".join([f"<tr><td>{ab.name}</td><td>{ab.fluorochrome}</td><td>{ab.target}</td><td>{ab.concentration} μg/mL</td></tr>" for ab in st.session_state.antibodies.values()])}
            </table>
            
            <h2>实验步骤摘要</h2>
            <ol>
                <li>样本准备与分装</li>
                <li>Fc阻断与表面染色 (30min, 4°C)</li>
                <li>固定与破膜 (20min, RT)</li>
                <li>胞内染色 (45min, 4°C)</li>
                <li>洗涤重悬上机</li>
            </ol>
            
            <h2>注意事项</h2>
            <ul>
                <li>所有操作避光进行</li>
                <li>离心条件: 300g, 4°C, 5分钟</li>
                <li>抗体现配现用</li>
                <li>设置正确的对照</li>
            </ul>
        </body>
        </html>
        """
        
        st.download_button(
            label="📘 下载实验方案 (HTML)",
            data=html_protocol,
            file_name=f"{st.session_state.current_project}_protocol.html",
            mime="text/html",
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
        st.markdown(f"**管子配置:** {len(st.session_state.tubes)}种")
        
        st.markdown("---")
        
        # 导航菜单
        st.markdown("### 📋 导航菜单")
        
        page_options = {
            "🏠 仪表盘": "dashboard",
            "🧪 抗体管理": "antibodies",
            "🧫 管子设计": "tubes",
            "🔢 染色矩阵": "matrix",
            "🧪 母液计算": "mastermix",
            "📋 实验计划": "planner",
            "📖 实验方案": "protocol"
        }
        
        selected_page = st.radio(
            "选择页面",
            options=list(page_options.keys()),
            label_visibility="collapsed"
        )
        
        st.markdown("---")
        
        # 快速操作
        st.markdown("### ⚡ 快速操作")
        
        col_q1, col_q2 = st.columns(2)
        with col_q1:
            if st.button("🔄 重置", use_container_width=True, type="secondary"):
                for key in ['antibodies', 'tubes']:
                    if key in st.session_state:
                        st.session_state[key] = {}
                st.rerun()
        
        with col_q2:
            if st.button("💾 保存项目", use_container_width=True):
                project_data = {
                    'project_name': st.session_state.current_project,
                    'antibodies': {name: asdict(ab) for name, ab in st.session_state.antibodies.items()},
                    'tubes': {name: asdict(tube) for name, tube in st.session_state.tubes.items()},
                    'volumes': st.session_state.volumes
                }
                
                json_str = json.dumps(project_data, ensure_ascii=False, indent=2, default=str)
                st.download_button(
                    label="下载JSON",
                    data=json_str,
                    file_name=f"{st.session_state.current_project}.json",
                    mime="application/json",
                    use_container_width=True
                )
        
        # 加载项目
        uploaded_file = st.file_uploader("📂 加载项目文件", type=['json'])
        if uploaded_file is not None:
            try:
                project_data = json.load(uploaded_file)
                
                # 加载抗体
                st.session_state.antibodies.clear()
                for name, ab_dict in project_data.get('antibodies', {}).items():
                    # 转换类型字符串为枚举
                    type_map = {t.value: t for t in AntibodyType}
                    ab_dict['type'] = type_map[ab_dict['type']]
                    st.session_state.antibodies[name] = Antibody(**ab_dict)
                
                # 加载管子
                st.session_state.tubes.clear()
                for name, tube_dict in project_data.get('tubes', {}).items():
                    st.session_state.tubes[name] = TubeConfiguration(**tube_dict)
                
                st.session_state.volumes = project_data.get('volumes', st.session_state.volumes)
                st.session_state.current_project = project_data.get('project_name', st.session_state.current_project)
                
                st.success("✅ 项目加载成功！")
                st.rerun()
                
            except Exception as e:
                st.error(f"加载失败: {e}")
        
        st.markdown("---")
        st.markdown("### ℹ️ 使用说明")
        st.markdown("""
        1. 从"抗体管理"开始添加抗体
        2. 在"管子设计"中配置实验管子
        3. 查看"染色矩阵"确认配置
        4. 使用"母液计算"获取配方
        5. 生成"实验计划"和"实验方案"
        """)
        
        return page_options[selected_page]

def render_dashboard():
    """仪表盘页面"""
    st.markdown('<div class="main-header">🔬 流式染色矩阵设计器</div>', unsafe_allow_html=True)
    
    # 项目概览
    col_dash1, col_dash2, col_dash3, col_dash4 = st.columns(4)
    
    with col_dash1:
        st.metric("抗体数量", len(st.session_state.antibodies))
    with col_dash2:
        st.metric("管子配置", len(st.session_state.tubes))
    with col_dash3:
        # 计算表面抗体数量
        surface_count = sum(1 for ab in st.session_state.antibodies.values() 
                          if ab.type == AntibodyType.SURFACE)
        st.metric("表面抗体", surface_count)
    with col_dash4:
        # 计算胞内抗体数量
        intracel_count = sum(1 for ab in st.session_state.antibodies.values() 
                           if ab.type == AntibodyType.INTRACELLULAR)
        st.metric("胞内抗体", intracel_count)
    
    # 快速开始指南
    st.markdown("### 🚀 快速开始")
    
    guide_col1, guide_col2 = st.columns(2)
    
    with guide_col1:
        st.markdown("""
        #### 新手模式
        1. 点击"加载标准实验模板"
        2. 查看生成的染色矩阵
        3. 使用母液计算器获取配方
        4. 导出实验计划
        
        **适合**: CD45 + α-SMA实验
        """)
        
        if st.button("📥 加载标准模板", use_container_width=True):
            designer = FlowMatrixDesigner()
            designer.load_standard_antibodies()
            designer.load_standard_tubes()
            st.success("✅ 标准模板加载成功！")
            st.rerun()
    
    with guide_col2:
        st.markdown("""
        #### 自定义模式
        1. 在"抗体管理"中添加抗体
        2. 在"管子设计"中配置实验管
        3. 查看和调整染色矩阵
        4. 生成完整的实验方案
        
        **适合**: 自定义多色panel
        """)
        
        if st.button("🆕 开始自定义设计", use_container_width=True):
            st.switch_page("抗体管理")
    
    # 最近活动
    if st.session_state.antibodies or st.session_state.tubes:
        st.markdown("### 📊 当前项目概览")
        
        col_overview1, col_overview2 = st.columns(2)
        
        with col_overview1:
            if st.session_state.antibodies:
                st.markdown("#### 抗体列表")
                for ab in list(st.session_state.antibodies.values())[:5]:  # 显示前5个
                    st.markdown(f"- **{ab.short_name}**: {ab.target} ({ab.fluorochrome})")
                if len(st.session_state.antibodies) > 5:
                    st.caption(f"... 还有 {len(st.session_state.antibodies) - 5} 个抗体")
        
        with col_overview2:
            if st.session_state.tubes:
                st.markdown("#### 管子配置")
                for tube in list(st.session_state.tubes.values())[:5]:  # 显示前5个
                    st.markdown(f"- **{tube.name}**: {len(tube.antibodies)}种抗体")
                if len(st.session_state.tubes) > 5:
                    st.caption(f"... 还有 {len(st.session_state.tubes) - 5} 种管子")
    
    # 功能卡片
    st.markdown("### 🛠️ 功能模块")
    
    col_func1, col_func2, col_func3 = st.columns(3)
    
    with col_func1:
        st.markdown("""
        <div class="antibody-card">
            <h4>🧪 抗体管理</h4>
            <p>添加、编辑和管理抗体信息</p>
        </div>
        """, unsafe_allow_html=True)
        
        if st.button("前往抗体管理", key="goto_ab", use_container_width=True):
            st.switch_page("antibodies")
    
    with col_func2:
        st.markdown("""
        <div class="tube-card">
            <h4>🧫 管子设计</h4>
            <p>配置实验管和对照管</p>
        </div>
        """, unsafe_allow_html=True)
        
        if st.button("前往管子设计", key="goto_tube", use_container_width=True):
            st.switch_page("tubes")
    
    with col_func3:
        st.markdown("""
        <div style="
            background: linear-gradient(135deg, #10B98120 0%, #10B98140 100%);
            border: 2px solid #10B981;
            padding: 1rem;
            border-radius: 12px;
            margin: 0.5rem 0;
        ">
            <h4>🔢 染色矩阵</h4>
            <p>可视化染色方案矩阵</p>
        </div>
        """, unsafe_allow_html=True)
        
        if st.button("前往染色矩阵", key="goto_matrix", use_container_width=True):
            st.switch_page("matrix")

def main():
    """主函数"""
    # 初始化session state
    init_session()
    
    # 渲染侧边栏并获取当前页面
    current_page = render_sidebar()
    
    # 根据选择渲染页面
    if current_page == "dashboard":
        render_dashboard()
    elif current_page == "antibodies":
        render_antibody_management()
    elif current_page == "tubes":
        render_tube_design()
    elif current_page == "matrix":
        render_matrix()
    elif current_page == "mastermix":
        render_mastermix_calculator()
    elif current_page == "planner":
        render_experiment_planner()
    elif current_page == "protocol":
        render_protocol()

if __name__ == "__main__":
    main()
