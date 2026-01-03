# -*- coding: utf-8 -*-
"""
Created on Sat Jan  3 17:38:30 2026

@author: sjc
"""

#!/usr/bin/env python3
"""
Flow Cytometry Staining Matrix Designer - 交互式命令行工具
"""

import sys
import os
import json
from typing import Dict, List, Optional
from dataclasses import dataclass, field, asdict
from enum import Enum
import pandas as pd
from datetime import datetime

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
            self.short_name = self.name.split()[-1] if ' ' in self.name else self.name[:10]

@dataclass
class TubeConfiguration:
    """管子配置"""
    name: str
    description: str
    antibodies: List[str] = field(default_factory=list)
    needs_fixation: bool = False
    is_control: bool = False
    control_type: str = ""  # FMO, Isotype, Single, Blank
    
    def add_antibody(self, antibody_name: str):
        if antibody_name not in self.antibodies:
            self.antibodies.append(antibody_name)
    
    def remove_antibody(self, antibody_name: str):
        if antibody_name in self.antibodies:
            self.antibodies.remove(antibody_name)

class FlowMatrixDesigner:
    """流式矩阵设计器"""
    
    def __init__(self):
        self.antibodies: Dict[str, Antibody] = {}
        self.tubes: Dict[str, TubeConfiguration] = {}
        self.current_project: str = f"Flow_Project_{datetime.now().strftime('%Y%m%d_%H%M')}"
        self.volumes = {
            'per_tube': 100.0,
            'intracellular_per_tube': 50.0,
            'cell_count': 1.0,  # 10⁶ cells
            'extra_tubes': 2
        }
        
    def clear_screen(self):
        """清屏"""
        os.system('cls' if os.name == 'nt' else 'clear')
    
    def print_header(self, title: str):
        """打印标题"""
        self.clear_screen()
        print("=" * 60)
        print(f"流式染色矩阵设计器".center(60))
        print(f"{title}".center(60))
        print("=" * 60)
        print()
    
    def wait_for_continue(self):
        """等待用户继续"""
        input("\n按 Enter 继续...")
    
    def show_main_menu(self):
        """显示主菜单"""
        self.print_header("主菜单")
        print("请选择操作:")
        print(" 1. 📋 管理抗体库")
        print(" 2. 🧪 设计实验管子")
        print(" 3. 🔢 生成染色矩阵")
        print(" 4. 🧪 计算母液配方")
        print(" 5. 📊 生成实验计划")
        print(" 6. 💾 保存项目")
        print(" 7. 📂 加载项目")
        print(" 8. 🖨️ 打印实验方案")
        print(" 9. 🚪 退出")
        print()
        
        choice = input("请输入选择 (1-9): ").strip()
        return choice
    
    def manage_antibodies(self):
        """管理抗体库"""
        while True:
            self.print_header("抗体库管理")
            print(f"当前有 {len(self.antibodies)} 种抗体")
            print()
            print("选择操作:")
            print(" 1. 添加新抗体")
            print(" 2. 查看抗体列表")
            print(" 3. 编辑抗体")
            print(" 4. 删除抗体")
            print(" 5. 导入标准抗体（CD45+α-SMA实验）")
            print(" 6. 返回主菜单")
            print()
            
            choice = input("请选择 (1-6): ").strip()
            
            if choice == '1':
                self.add_antibody()
            elif choice == '2':
                self.list_antibodies()
            elif choice == '3':
                self.edit_antibody()
            elif choice == '4':
                self.delete_antibody()
            elif choice == '5':
                self.load_standard_antibodies()
            elif choice == '6':
                break
            else:
                print("无效选择！")
                self.wait_for_continue()
    
    def add_antibody(self):
        """添加抗体"""
        self.print_header("添加新抗体")
        
        print("请输入抗体信息（按Enter跳过可选字段）:")
        print()
        
        name = input("抗体全名: ").strip()
        if not name:
            print("抗体名称不能为空！")
            self.wait_for_continue()
            return
        
        short_name = input(f"简称 [{name.split()[-1] if ' ' in name else name[:10]}]: ").strip()
        fluorochrome = input("荧光染料 (如 AF647, BB515): ").strip()
        target = input("靶标 (如 CD45, α-SMA): ").strip()
        clone = input("克隆号: ").strip()
        
        try:
            concentration = float(input("浓度 (μg/mL): ").strip() or "0")
            recommended_use = float(input("推荐用量 (μg/10⁶ cells): ").strip() or "0")
        except ValueError:
            print("浓度和用量必须是数字！")
            self.wait_for_continue()
            return
        
        print("\n抗体类型:")
        print(" 1. 表面抗体")
        print(" 2. 胞内抗体")
        print(" 3. 死活染料")
        print(" 4. Fc阻断剂")
        print(" 5. 其他")
        
        type_choice = input("选择类型 (1-5): ").strip()
        type_map = {
            '1': AntibodyType.SURFACE,
            '2': AntibodyType.INTRACELLULAR,
            '3': AntibodyType.VIABILITY,
            '4': AntibodyType.FC_BLOCK,
            '5': AntibodyType.OTHER
        }
        
        catalog = input("货号: ").strip()
        lot = input("批号: ").strip()
        notes = input("备注: ").strip()
        
        antibody = Antibody(
            name=name,
            short_name=short_name or name.split()[-1] if ' ' in name else name[:10],
            fluorochrome=fluorochrome,
            target=target,
            clone=clone,
            concentration=concentration,
            recommended_use=recommended_use,
            type=type_map.get(type_choice, AntibodyType.SURFACE),
            catalog_number=catalog,
            lot_number=lot,
            notes=notes
        )
        
        self.antibodies[antibody.name] = antibody
        print(f"\n✅ 已成功添加抗体: {antibody.name}")
        self.wait_for_continue()
    
    def list_antibodies(self):
        """列出所有抗体"""
        self.print_header("抗体列表")
        
        if not self.antibodies:
            print("尚无抗体信息")
            self.wait_for_continue()
            return
        
        for i, (name, ab) in enumerate(self.antibodies.items(), 1):
            print(f"{i:2d}. {name}")
            print(f"     靶标: {ab.target}, 荧光: {ab.fluorochrome}, 类型: {ab.type.value}")
            print(f"     浓度: {ab.concentration} μg/mL, 用量: {ab.recommended_use} μg/10⁶ cells")
            print(f"     货号: {ab.catalog_number}, 批号: {ab.lot_number}")
            if ab.notes:
                print(f"     备注: {ab.notes}")
            print()
        
        self.wait_for_continue()
    
    def edit_antibody(self):
        """编辑抗体"""
        if not self.antibodies:
            print("尚无抗体可编辑")
            self.wait_for_continue()
            return
        
        self.list_antibodies()
        try:
            choice = int(input("\n选择要编辑的抗体编号: ").strip())
            ab_list = list(self.antibodies.keys())
            if 1 <= choice <= len(ab_list):
                ab_name = ab_list[choice-1]
                print(f"\n编辑抗体: {ab_name}")
                # 这里可以添加详细的编辑逻辑
                print("编辑功能开发中...")
            else:
                print("无效编号！")
        except ValueError:
            print("请输入数字！")
        
        self.wait_for_continue()
    
    def delete_antibody(self):
        """删除抗体"""
        if not self.antibodies:
            print("尚无抗体可删除")
            self.wait_for_continue()
            return
        
        self.list_antibodies()
        try:
            choice = int(input("\n选择要删除的抗体编号: ").strip())
            ab_list = list(self.antibodies.keys())
            if 1 <= choice <= len(ab_list):
                ab_name = ab_list[choice-1]
                confirm = input(f"确认删除抗体 '{ab_name}'? (y/N): ").strip().lower()
                if confirm == 'y':
                    del self.antibodies[ab_name]
                    print(f"✅ 已删除抗体: {ab_name}")
            else:
                print("无效编号！")
        except ValueError:
            print("请输入数字！")
        
        self.wait_for_continue()
    
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
                notes="Fc受体阻断剂"
            ),
            Antibody(
                name="Live/Dye eF780",
                short_name="LiveDye",
                fluorochrome="eF780",
                target="Viability",
                clone="N/A",
                concentration=1000.0,
                recommended_use=0.5,  # 1:1000稀释
                type=AntibodyType.VIABILITY,
                catalog_number="65-0865-14",
                lot_number="123456",
                notes="死活染料，建议1:1000稀释使用"
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
                target="α-Smooth Muscle Actin",
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
            self.antibodies[ab.name] = ab
        
        print("✅ 已加载标准抗体（CD45+α-SMA实验）")
        self.wait_for_continue()
    
    def design_tubes(self):
        """设计实验管子"""
        while True:
            self.print_header("实验管子设计")
            print(f"当前有 {len(self.tubes)} 种管子配置")
            print()
            print("选择操作:")
            print(" 1. 创建新管子配置")
            print(" 2. 查看管子列表")
            print(" 3. 编辑管子配置")
            print(" 4. 删除管子")
            print(" 5. 导入标准配置（基于您的矩阵）")
            print(" 6. 返回主菜单")
            print()
            
            choice = input("请选择 (1-6): ").strip()
            
            if choice == '1':
                self.create_tube()
            elif choice == '2':
                self.list_tubes()
            elif choice == '3':
                self.edit_tube()
            elif choice == '4':
                self.delete_tube()
            elif choice == '5':
                self.load_standard_tubes()
            elif choice == '6':
                break
            else:
                print("无效选择！")
                self.wait_for_continue()
    
    def create_tube(self):
        """创建管子配置"""
        self.print_header("创建新管子配置")
        
        name = input("管子名称 (如: Full_Stain, FMO_αSMA): ").strip()
        if not name:
            print("名称不能为空！")
            self.wait_for_continue()
            return
        
        description = input("描述: ").strip()
        
        print("\n是否为对照管?")
        is_control = input("是对照管吗? (y/N): ").strip().lower() == 'y'
        
        control_type = ""
        if is_control:
            print("\n对照类型:")
            print(" 1. FMO (荧光减一)")
            print(" 2. 同型对照")
            print(" 3. 单阳对照")
            print(" 4. 空白对照")
            ct_choice = input("选择类型 (1-4): ").strip()
            control_map = {'1': 'FMO', '2': 'Isotype', '3': 'Single', '4': 'Blank'}
            control_type = control_map.get(ct_choice, '')
        
        needs_fixation = input("需要固定破膜吗? (y/N): ").strip().lower() == 'y'
        
        tube = TubeConfiguration(
            name=name,
            description=description,
            needs_fixation=needs_fixation,
            is_control=is_control,
            control_type=control_type
        )
        
        # 添加抗体
        if self.antibodies:
            print("\n选择要添加到管子的抗体:")
            ab_list = list(self.antibodies.keys())
            for i, ab_name in enumerate(ab_list, 1):
                print(f" {i}. {ab_name}")
            
            while True:
                choice = input("\n输入抗体编号 (输入0完成): ").strip()
                if choice == '0':
                    break
                try:
                    idx = int(choice) - 1
                    if 0 <= idx < len(ab_list):
                        tube.add_antibody(ab_list[idx])
                        print(f"✅ 已添加: {ab_list[idx]}")
                    else:
                        print("无效编号！")
                except ValueError:
                    print("请输入数字！")
        
        self.tubes[name] = tube
        print(f"\n✅ 已创建管子: {name}")
        self.wait_for_continue()
    
    def list_tubes(self):
        """列出所有管子"""
        self.print_header("管子配置列表")
        
        if not self.tubes:
            print("尚无管子配置")
            self.wait_for_continue()
            return
        
        for name, tube in self.tubes.items():
            print(f"🔬 {name}")
            print(f"   描述: {tube.description}")
            print(f"   抗体: {', '.join(tube.antibodies) if tube.antibodies else '无'}")
            print(f"   固定破膜: {'是' if tube.needs_fixation else '否'}")
            if tube.is_control:
                print(f"   对照类型: {tube.control_type}")
            print()
        
        self.wait_for_continue()
    
    def load_standard_tubes(self):
        """加载标准管子配置（基于您的矩阵）"""
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
        
        self.tubes = standard_tubes
        print("✅ 已加载标准管子配置")
        self.wait_for_continue()
    
    def generate_matrix(self):
        """生成染色矩阵"""
        self.print_header("染色矩阵")
        
        if not self.tubes or not self.antibodies:
            print("需要先配置抗体和管子！")
            self.wait_for_continue()
            return
        
        # 创建矩阵数据
        matrix_data = []
        ab_names = list(self.antibodies.keys())
        
        for tube_name, tube in self.tubes.items():
            row = {"管子名称": tube_name, "描述": tube.description}
            for ab_name in ab_names:
                row[ab_name] = "✓" if ab_name in tube.antibodies else ""
            row["固定破膜"] = "是" if tube.needs_fixation else "否"
            row["对照类型"] = tube.control_type if tube.is_control else "实验管"
            matrix_data.append(row)
        
        # 创建DataFrame
        df = pd.DataFrame(matrix_data)
        columns = ["管子名称", "描述"] + ab_names + ["固定破膜", "对照类型"]
        df = df[columns]
        
        # 显示矩阵
        pd.set_option('display.max_columns', None)
        pd.set_option('display.width', None)
        pd.set_option('display.max_colwidth', 30)
        
        print(df.to_string(index=False))
        print()
        
        # 导出选项
        export = input("是否导出为CSV文件? (y/N): ").strip().lower()
        if export == 'y':
            filename = f"{self.current_project}_matrix.csv"
            df.to_csv(filename, index=False, encoding='utf-8-sig')
            print(f"✅ 矩阵已导出到: {filename}")
        
        self.wait_for_continue()
    
    def calculate_mastermix(self):
        """计算母液配方"""
        self.print_header("母液配方计算")
        
        if not self.antibodies:
            print("需要先配置抗体！")
            self.wait_for_continue()
            return
        
        print("实验参数设置:")
        try:
            cell_count = float(input(f"细胞数 (×10⁶) [{self.volumes['cell_count']}]: ").strip() or str(self.volumes['cell_count']))
            per_tube = float(input(f"每管体积 (μL) [{self.volumes['per_tube']}]: ").strip() or str(self.volumes['per_tube']))
            extra_tubes = int(input(f"安全余量管数 [{self.volumes['extra_tubes']}]: ").strip() or str(self.volumes['extra_tubes']))
            
            self.volumes.update({
                'cell_count': cell_count,
                'per_tube': per_tube,
                'extra_tubes': extra_tubes
            })
        except ValueError:
            print("输入无效，使用默认值")
        
        print("\n" + "="*60)
        print("母液配方计算结果")
        print("="*60)
        
        # 统计需要不同混合液的管子
        surface_tubes = []
        intracellular_tubes = []
        
        for tube_name, tube in self.tubes.items():
            if tube.needs_fixation:
                intracellular_tubes.append(tube_name)
            elif tube.antibodies:  # 有抗体且不需要破膜
                surface_tubes.append(tube_name)
        
        # 1. 表面染色母液
        if surface_tubes:
            print(f"\n🔬 表面染色母液 (用于 {len(surface_tubes)} 管)")
            print(f"   适用管子: {', '.join(surface_tubes)}")
            
            total_tubes = len(surface_tubes) + extra_tubes
            total_volume = per_tube * total_tubes
            
            print(f"   总体积: {total_volume} μL (每管{per_tube}μL × {total_tubes}管)")
            print(f"   配制步骤:")
            print(f"     1. 取洁净EP管")
            print(f"     2. 加入以下抗体:")
            
            for tube_name in surface_tubes:
                tube = self.tubes[tube_name]
                for ab_name in tube.antibodies:
                    if ab_name in self.antibodies:
                        ab = self.antibodies[ab_name]
                        if ab.type in [AntibodyType.SURFACE, AntibodyType.VIABILITY, AntibodyType.FC_BLOCK]:
                            per_tube_vol = (ab.recommended_use * cell_count) / ab.concentration
                            total_vol = per_tube_vol * total_tubes
                            print(f"       - {ab.name}")
                            print(f"         每管: {per_tube_vol:.2f} μL, 总计: {total_vol:.2f} μL")
            
            print(f"     3. 用流式染色缓冲液补至 {total_volume} μL")
            print(f"     4. 涡旋混匀，避光保存")
        
        # 2. 胞内染色工作液
        if intracellular_tubes:
            print(f"\n🧫 胞内染色工作液 (用于 {len(intracellular_tubes)} 管)")
            print(f"   适用管子: {', '.join(intracellular_tubes)}")
            
            intracel_volume = self.volumes['intracellular_per_tube']
            total_tubes = len(intracellular_tubes) + extra_tubes
            total_volume = intracel_volume * total_tubes
            
            print(f"   总体积: {total_volume} μL (每管{intracel_volume}μL × {total_tubes}管)")
            print(f"   用1X破膜缓冲液配制:")
            
            for tube_name in intracellular_tubes:
                tube = self.tubes[tube_name]
                for ab_name in tube.antibodies:
                    if ab_name in self.antibodies:
                        ab = self.antibodies[ab_name]
                        if ab.type == AntibodyType.INTRACELLULAR:
                            per_tube_vol = (ab.recommended_use * cell_count) / ab.concentration
                            total_vol = per_tube_vol * total_tubes
                            print(f"       - {ab.name}")
                            print(f"         每管: {per_tube_vol:.2f} μL, 总计: {total_vol:.2f} μL")
            
            print(f"   用1X破膜缓冲液补至 {total_volume} μL")
            print(f"   涡旋混匀，4°C避光保存")
        
        print()
        self.wait_for_continue()
    
    def generate_experiment_plan(self):
        """生成实验计划"""
        self.print_header("实验计划生成")
        
        if not self.tubes:
            print("需要先配置管子！")
            self.wait_for_continue()
            return
        
        print("设置实验组:")
        groups_input = input("输入实验组名称（用逗号分隔，如：Control,Model,Treatment）: ").strip()
        groups = [g.strip() for g in groups_input.split(',')] if groups_input else ["Control", "Model", "Treatment"]
        
        try:
            replicates = int(input("每组重复数 (默认3): ").strip() or "3")
        except ValueError:
            replicates = 3
        
        print("\n" + "="*60)
        print("实验计划")
        print("="*60)
        
        plan_data = []
        sample_id = 1
        
        for group in groups:
            for rep in range(1, replicates + 1):
                for tube_name, tube in self.tubes.items():
                    plan_data.append({
                        "样品ID": f"S{sample_id:03d}",
                        "实验组": group,
                        "重复": rep,
                        "管子类型": tube_name,
                        "描述": tube.description,
                        "抗体": ", ".join(tube.antibodies),
                        "固定破膜": "是" if tube.needs_fixation else "否",
                        "上机顺序": sample_id
                    })
                    sample_id += 1
        
        df = pd.DataFrame(plan_data)
        print(df.to_string(index=False))
        print(f"\n总计: {len(plan_data)} 个样品")
        
        # 导出
        export = input("\n是否导出实验计划? (y/N): ").strip().lower()
        if export == 'y':
            filename = f"{self.current_project}_plan.csv"
            df.to_csv(filename, index=False, encoding='utf-8-sig')
            print(f"✅ 实验计划已导出到: {filename}")
            
            # 同时生成上机列表
            run_order = df[["样品ID", "实验组", "管子类型", "上机顺序"]].copy()
            run_order_filename = f"{self.current_project}_run_order.csv"
            run_order.to_csv(run_order_filename, index=False)
            print(f"✅ 上机顺序表已导出到: {run_order_filename}")
        
        self.wait_for_continue()
    
    def save_project(self):
        """保存项目"""
        self.print_header("保存项目")
        
        filename = input(f"输入保存文件名 [{self.current_project}.json]: ").strip()
        if not filename:
            filename = f"{self.current_project}.json"
        elif not filename.endswith('.json'):
            filename += '.json'
        
        project_data = {
            'project_name': self.current_project,
            'antibodies': {name: asdict(ab) for name, ab in self.antibodies.items()},
            'tubes': {name: asdict(tube) for name, tube in self.tubes.items()},
            'volumes': self.volumes,
            'save_time': datetime.now().isoformat()
        }
        
        try:
            with open(filename, 'w', encoding='utf-8') as f:
                json.dump(project_data, f, ensure_ascii=False, indent=2)
            print(f"✅ 项目已保存到: {filename}")
        except Exception as e:
            print(f"❌ 保存失败: {e}")
        
        self.wait_for_continue()
    
    def load_project(self):
        """加载项目"""
        self.print_header("加载项目")
        
        import glob
        json_files = glob.glob("*.json")
        if json_files:
            print("可用的项目文件:")
            for i, f in enumerate(json_files, 1):
                print(f" {i}. {f}")
            print()
        
        filename = input("输入要加载的文件名: ").strip()
        if not filename.endswith('.json'):
            filename += '.json'
        
        try:
            with open(filename, 'r', encoding='utf-8') as f:
                project_data = json.load(f)
            
            # 加载抗体
            self.antibodies.clear()
            for name, ab_dict in project_data.get('antibodies', {}).items():
                ab_dict['type'] = AntibodyType(ab_dict['type'])  # 转换回枚举
                self.antibodies[name] = Antibody(**ab_dict)
            
            # 加载管子
            self.tubes.clear()
            for name, tube_dict in project_data.get('tubes', {}).items():
                self.tubes[name] = TubeConfiguration(**tube_dict)
            
            self.volumes = project_data.get('volumes', self.volumes)
            self.current_project = project_data.get('project_name', self.current_project)
            
            print(f"✅ 项目已加载: {self.current_project}")
            print(f"   抗体: {len(self.antibodies)} 种")
            print(f"   管子配置: {len(self.tubes)} 种")
            
        except FileNotFoundError:
            print(f"❌ 文件不存在: {filename}")
        except Exception as e:
            print(f"❌ 加载失败: {e}")
        
        self.wait_for_continue()
    
    def print_protocol(self):
        """打印实验方案"""
        self.print_header("实验方案")
        
        protocol = f"""
        流式染色实验方案
        项目: {self.current_project}
        生成时间: {datetime.now().strftime('%Y-%m-%d %H:%M')}
        
        =============== 实验信息 ===============
        抗体种类: {len(self.antibodies)}
        管子配置: {len(self.tubes)}
        
        =============== 染色步骤 ===============
        1. 样本准备
           - 制备单细胞悬液，调整浓度至 1×10⁷ cells/mL
           - 按实验计划分装细胞到标记好的流式管中
           - 每管分装 1×10⁶ cells (100 μL)
        
        2. Fc受体阻断与表面染色
           a. 配制表面染色母液
           b. 向对应管子中加入 {self.volumes['per_tube']} μL 母液
           c. 4°C避光孵育30分钟
           d. 加入1 mL预冷染色缓冲液，300g 4°C离心5分钟
           e. 弃上清，重复洗涤一次
        
        3. 固定与破膜 (仅需胞内染色的管子)
           a. 加入100 μL固定液
           b. 室温避光孵育20分钟
           c. 加入1 mL 1X破膜缓冲液，300g 4°C离心5分钟
           d. 弃上清，重复洗涤一次
        
        4. 胞内染色
           a. 用破膜缓冲液配制胞内抗体工作液
           b. 向对应管子中加入50 μL工作液
           c. 4°C避光孵育45分钟
           d. 用破膜缓冲液洗涤2次
        
        5. 上机检测
           a. 所有管子用300 μL流式染色缓冲液重悬
           b. 过35μm细胞筛网
           c. 按上机顺序表进行检测
        
        =============== 注意事项 ===============
        1. 所有操作避光进行
        2. 离心条件: 300g, 4°C, 5分钟
        3. 抗体母液现配现用
        4. 上机前充分混匀样本
        5. 设置正确的补偿对照
        """
        
        print(protocol)
        
        export = input("\n是否导出实验方案? (y/N): ").strip().lower()
        if export == 'y':
            filename = f"{self.current_project}_protocol.txt"
            with open(filename, 'w', encoding='utf-8') as f:
                f.write(protocol)
            print(f"✅ 实验方案已导出到: {filename}")
        
        self.wait_for_continue()
    
    def run(self):
        """运行主程序"""
        print("正在启动流式染色矩阵设计器...")
        
        # 检查是否需要加载标准配置
        if not self.antibodies:
            load_std = input("是否加载标准抗体和管子配置? (Y/n): ").strip().lower()
            if load_std in ['', 'y', 'yes']:
                self.load_standard_antibodies()
                self.load_standard_tubes()
        
        while True:
            try:
                choice = self.show_main_menu()
                
                if choice == '1':
                    self.manage_antibodies()
                elif choice == '2':
                    self.design_tubes()
                elif choice == '3':
                    self.generate_matrix()
                elif choice == '4':
                    self.calculate_mastermix()
                elif choice == '5':
                    self.generate_experiment_plan()
                elif choice == '6':
                    self.save_project()
                elif choice == '7':
                    self.load_project()
                elif choice == '8':
                    self.print_protocol()
                elif choice == '9':
                    print("\n感谢使用流式染色矩阵设计器！")
                    print("再见！👋")
                    break
                else:
                    print("无效选择，请重新输入！")
                    self.wait_for_continue()
                    
            except KeyboardInterrupt:
                print("\n\n程序被中断")
                save = input("是否保存当前项目? (y/N): ").strip().lower()
                if save == 'y':
                    self.save_project()
                break
            except Exception as e:
                print(f"\n❌ 发生错误: {e}")
                self.wait_for_continue()

def main():
    """主函数"""
    try:
        designer = FlowMatrixDesigner()
        designer.run()
    except Exception as e:
        print(f"程序运行出错: {e}")
        input("按 Enter 退出...")

if __name__ == "__main__":
    main()