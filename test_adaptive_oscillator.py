#!/usr/bin/env python3
"""
自适应振荡器测试文件
用于验证2阶自适应振荡器的相位估计性能
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

# 设置中文字体
plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'SimHei', 'WenQuanYi Micro Hei']
plt.rcParams['axes.unicode_minus'] = False

class AdaptiveOscillator:
    """2阶自适应振荡器类 - 基于MATLAB实现"""
    
    def __init__(self, initial_frequency=1.0, dt=0.02):
        self.dt = dt
        
        # 振荡器阶数
        self.order = 2
        
        # 振荡器参数 [幅值, 频率, 相位]
        self.para = np.zeros((self.order, 3))
        self.para_dot = np.zeros((self.order, 3))
        
        # 初始化参数
        self.para[0, 0] = 0.0
        self.para[0, 2] = np.pi / 2
        self.para[1, 0] = 0.4
        self.para[1, 1] = 2 * np.pi * initial_frequency
        self.para[1, 2] = 0.0
        
        # 自适应参数
        self.init_frequency = initial_frequency
        T = 1.0 / initial_frequency
        T_aefa = 2 * T
        T_amig = 2 * T
        
        self.yita = 2.0 / T_aefa
        self.vaomig = 20.0 / (T_amig ** 2)
        self.vfai = np.sqrt(24.2 * self.vaomig)
        
        # 振荡器输出
        self.y = 0.0
        
        # 频率限制
        self.omega_min = 2 * np.pi * 0.3
        self.omega_max = 2 * np.pi * 3.0
        
        # 相位输出
        self.current_phase = 0.0
        
    def _wrap_to_2pi(self, angle):
        """将角度标准化到 [0, 2π) 范围"""
        return angle % (2 * np.pi)
    
    def initialize(self, q, dq, frequency):
        """步态开始时初始化振荡器参数"""
        T = 1.0 / frequency if frequency > 0 else 1.0
        T_aefa = 2 * T
        T_amig = 2 * T
        
        self.yita = 2.0 / T_aefa
        self.vaomig = 20.0 / (T_amig ** 2)
        self.vfai = np.sqrt(24.2 * self.vaomig)
        
        self.para[0, 0] = 0.0
        self.para[0, 2] = np.pi / 2
        self.para[1, 0] = 0.4
        self.para[1, 1] = 2 * np.pi * frequency
        
        omega = self.para[1, 1]
        self.para[1, 2] = np.arctan2(omega * q, dq) if dq != 0 else 0.0
        
        sin_phi = np.sin(self.para[1, 2])
        if abs(sin_phi) > 0.01:
            self.para[1, 0] = max(q / sin_phi, 0.3)
        else:
            self.para[1, 0] = 0.4
        
        print(f"振荡器初始化 - 频率: {frequency:.2f} Hz, 相位: {self.para[1, 2]:.3f} rad, 幅值: {self.para[1, 0]:.3f}")
    
    def update(self, q, estimated_frequency=None):
        """更新振荡器状态"""
        # 计算振荡器输出
        q_esti = np.zeros(self.order)
        for j in range(self.order):
            q_esti[j] = self.para[j, 0] * np.sin(self.para[j, 2])
        
        self.y = np.sum(q_esti)
        F = q - self.y  # 误差信号
        
        # 更新参数变化率
        self.para_dot[0, 0] = self.yita * F
        
        for j in range(1, self.order):
            self.para_dot[j, 0] = self.yita * F * np.sin(self.para[j, 2])
        
        self.para_dot[1, 1] = self.vaomig * F * np.cos(self.para[1, 2])
        
        for j in range(1, self.order):
            self.para_dot[j, 2] = self.para[j, 1] + self.vfai * F * np.cos(self.para[j, 2])
            if self.para_dot[j, 2] < 0:
                self.para_dot[j, 2] = 0.0
        
        # 更新参数
        self.para[0, 0] = self.para[0, 0] + self.para_dot[0, 0] * self.dt
        self.para[0, 0] = 0.0
        
        for j in range(1, self.order):
            self.para[j, 0] = self.para[j, 0] + self.para_dot[j, 0] * self.dt
        
        self.para[1, 1] = self.para[1, 1] + self.para_dot[1, 1] * self.dt
        if self.para[1, 1] < 0:
            self.para[1, 1] = 0.0
        
        # 频率快速跟踪
        if estimated_frequency is not None:
            w_e = 2 * np.pi * estimated_frequency - self.para[1, 1]
            self.para[1, 1] = self.para[1, 1] + 5 * w_e * self.dt
        
        self.para[1, 1] = max(min(self.para[1, 1], self.omega_max), self.omega_min)
        
        for j in range(1, self.order):
            self.para[j, 2] = self.para[j, 2] + self.para_dot[j, 2] * self.dt
        
        phi = self._wrap_to_2pi(self.para[1, 2])
        self.para[1, 2] = phi
        self.current_phase = float(phi)
        
        return self.current_phase


def generate_gait_signal(t, frequency, amplitude=0.3, phase_shift=0):
    """
    生成模拟步态关节角度信号
    
    Args:
        t: 时间序列
        frequency: 步态频率 (Hz)
        amplitude: 信号幅值
        phase_shift: 相位偏移
    
    Returns:
        关节角度信号
    """
    omega = 2 * np.pi * frequency
    # 使用正弦波模拟关节角度
    angle = amplitude * np.sin(omega * t + phase_shift)
    return angle


def test_constant_frequency():
    """测试1: 恒定频率步态"""
    print("\n" + "="*60)
    print("测试1: 恒定频率步态 (1.0 Hz)")
    print("="*60)
    
    # 仿真参数
    dt = 0.02  # 采样时间 20ms (50Hz)
    duration = 10.0  # 仿真时长 10秒
    t = np.arange(0, duration, dt)
    n_samples = len(t)
    
    # 步态参数
    gait_frequency = 1.0  # Hz
    amplitude_left = 0.3
    amplitude_right = 0.3
    
    # 生成左右腿关节角度（相位差π）
    angle_left = generate_gait_signal(t, gait_frequency, amplitude_left, 0)
    angle_right = generate_gait_signal(t, gait_frequency, amplitude_right, np.pi)
    
    # 计算角度差（振荡器输入）
    angle_diff = angle_left - angle_right
    
    # 计算真实相位（用于对比）
    true_phase = np.mod(2 * np.pi * gait_frequency * t, 2 * np.pi)
    
    # 初始化振荡器
    oscillator = AdaptiveOscillator(initial_frequency=1.0, dt=dt)
    
    # 初始化时刻
    init_idx = 10  # 在第10个采样点初始化
    dq_init = (angle_diff[init_idx] - angle_diff[init_idx-1]) / dt
    oscillator.initialize(angle_diff[init_idx], dq_init, gait_frequency)
    
    # 存储结果
    estimated_phase = np.zeros(n_samples)
    estimated_frequency = np.zeros(n_samples)
    estimated_amplitude = np.zeros(n_samples)
    oscillator_output = np.zeros(n_samples)
    phase_error = np.zeros(n_samples)
    
    # 运行振荡器
    for i in range(n_samples):
        if i >= init_idx:
            estimated_phase[i] = oscillator.update(angle_diff[i], gait_frequency)
            estimated_frequency[i] = oscillator.para[1, 1] / (2 * np.pi)
            estimated_amplitude[i] = oscillator.para[1, 0]
            oscillator_output[i] = oscillator.y
            
            # 计算相位误差
            phase_error[i] = np.abs(np.mod(estimated_phase[i] - true_phase[i] + np.pi, 2*np.pi) - np.pi)
    
    # 绘图
    fig = plt.figure(figsize=(14, 10))
    gs = GridSpec(4, 2, figure=fig)
    
    # 1. 关节角度信号
    ax1 = fig.add_subplot(gs[0, :])
    ax1.plot(t, angle_left, 'b-', label='Left Leg', linewidth=1.5)
    ax1.plot(t, angle_right, 'r-', label='Right Leg', linewidth=1.5)
    ax1.plot(t, angle_diff, 'g-', label='Angle Diff (Input)', linewidth=1.5, alpha=0.7)
    ax1.axvline(x=t[init_idx], color='k', linestyle='--', label='Initialize', alpha=0.5)
    ax1.set_ylabel('Angle (rad)')
    ax1.set_title('Test 1: Constant Frequency Gait (1.0 Hz)')
    ax1.legend(loc='upper right')
    ax1.grid(True, alpha=0.3)
    
    # 2. 振荡器输出 vs 输入
    ax2 = fig.add_subplot(gs[1, :])
    ax2.plot(t, angle_diff, 'g-', label='Input Signal', linewidth=1.5)
    ax2.plot(t, oscillator_output, 'r--', label='Oscillator Output', linewidth=1.5)
    ax2.set_ylabel('Angle (rad)')
    ax2.set_title('Oscillator Tracking Performance')
    ax2.legend(loc='upper right')
    ax2.grid(True, alpha=0.3)
    
    # 3. 相位估计
    ax3 = fig.add_subplot(gs[2, 0])
    ax3.plot(t, true_phase, 'k-', label='True Phase', linewidth=2)
    ax3.plot(t, estimated_phase, 'b--', label='Estimated Phase', linewidth=1.5)
    ax3.set_ylabel('Phase (rad)')
    ax3.set_xlabel('Time (s)')
    ax3.set_title('Phase Estimation')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    
    # 4. 相位误差
    ax4 = fig.add_subplot(gs[2, 1])
    ax4.plot(t, np.rad2deg(phase_error), 'r-', linewidth=1.5)
    ax4.set_ylabel('Phase Error (deg)')
    ax4.set_xlabel('Time (s)')
    ax4.set_title('Phase Estimation Error')
    ax4.grid(True, alpha=0.3)
    
    # 5. 频率估计
    ax5 = fig.add_subplot(gs[3, 0])
    ax5.plot(t, estimated_frequency, 'b-', linewidth=1.5)
    ax5.axhline(y=gait_frequency, color='k', linestyle='--', label='True Frequency')
    ax5.set_ylabel('Frequency (Hz)')
    ax5.set_xlabel('Time (s)')
    ax5.set_title('Frequency Estimation')
    ax5.legend()
    ax5.grid(True, alpha=0.3)
    
    # 6. 幅值估计
    ax6 = fig.add_subplot(gs[3, 1])
    ax6.plot(t, estimated_amplitude, 'g-', linewidth=1.5)
    ax6.axhline(y=2*amplitude_left, color='k', linestyle='--', label='Expected Amplitude')
    ax6.set_ylabel('Amplitude')
    ax6.set_xlabel('Time (s)')
    ax6.set_title('Amplitude Estimation')
    ax6.legend()
    ax6.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('/home/zhang/gait_control_ws/test_constant_frequency.png', dpi=150)
    print(f"图像已保存: test_constant_frequency.png")
    
    # 统计信息
    converged_idx = int(2.0 / dt)  # 收敛后的数据（2秒后）
    if converged_idx < n_samples:
        mean_error = np.mean(phase_error[converged_idx:])
        std_error = np.std(phase_error[converged_idx:])
        max_error = np.max(phase_error[converged_idx:])
        
        print(f"\n收敛后统计 (2s后):")
        print(f"  平均相位误差: {np.rad2deg(mean_error):.3f}°")
        print(f"  相位误差标准差: {np.rad2deg(std_error):.3f}°")
        print(f"  最大相位误差: {np.rad2deg(max_error):.3f}°")


def test_frequency_change():
    """测试2: 变化频率步态"""
    print("\n" + "="*60)
    print("测试2: 变化频率步态 (0.8 Hz → 1.5 Hz)")
    print("="*60)
    
    dt = 0.02
    duration = 15.0
    t = np.arange(0, duration, dt)
    n_samples = len(t)
    
    # 生成变频信号
    angle_diff = np.zeros(n_samples)
    true_phase = np.zeros(n_samples)
    true_frequency = np.zeros(n_samples)
    
    change_time = 7.5  # 在7.5秒时改变频率
    freq1 = 0.8  # Hz
    freq2 = 1.5  # Hz
    amplitude = 0.3
    
    for i, ti in enumerate(t):
        if ti < change_time:
            freq = freq1
        else:
            freq = freq2
        
        true_frequency[i] = freq
        
        # 计算相位（保持连续性）
        if i == 0:
            true_phase[i] = 0
        else:
            true_phase[i] = true_phase[i-1] + 2 * np.pi * freq * dt
        
        # 生成信号
        angle_left = amplitude * np.sin(true_phase[i])
        angle_right = amplitude * np.sin(true_phase[i] + np.pi)
        angle_diff[i] = angle_left - angle_right
    
    # Wrap相位
    true_phase = np.mod(true_phase, 2 * np.pi)
    
    # 初始化振荡器
    oscillator = AdaptiveOscillator(initial_frequency=freq1, dt=dt)
    init_idx = 10
    dq_init = (angle_diff[init_idx] - angle_diff[init_idx-1]) / dt
    oscillator.initialize(angle_diff[init_idx], dq_init, freq1)
    
    # 运行振荡器
    estimated_phase = np.zeros(n_samples)
    estimated_frequency = np.zeros(n_samples)
    
    for i in range(n_samples):
        if i >= init_idx:
            # 使用真实频率进行快速跟踪
            estimated_phase[i] = oscillator.update(angle_diff[i], true_frequency[i])
            estimated_frequency[i] = oscillator.para[1, 1] / (2 * np.pi)
    
    # 绘图
    fig, axes = plt.subplots(3, 1, figsize=(12, 10))
    
    # 输入信号
    axes[0].plot(t, angle_diff, 'b-', linewidth=1.5)
    axes[0].axvline(x=change_time, color='r', linestyle='--', label='Frequency Change', linewidth=2)
    axes[0].axvline(x=t[init_idx], color='k', linestyle='--', label='Initialize', alpha=0.5)
    axes[0].set_ylabel('Angle Diff (rad)')
    axes[0].set_title('Test 2: Frequency Change (0.8 Hz → 1.5 Hz)')
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)
    
    # 频率估计
    axes[1].plot(t, true_frequency, 'k-', label='True Frequency', linewidth=2)
    axes[1].plot(t, estimated_frequency, 'b--', label='Estimated Frequency', linewidth=1.5)
    axes[1].axvline(x=change_time, color='r', linestyle='--', alpha=0.5)
    axes[1].set_ylabel('Frequency (Hz)')
    axes[1].set_title('Frequency Tracking')
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)
    
    # 相位估计
    axes[2].plot(t, true_phase, 'k-', label='True Phase', linewidth=2)
    axes[2].plot(t, estimated_phase, 'b--', label='Estimated Phase', linewidth=1.5)
    axes[2].axvline(x=change_time, color='r', linestyle='--', alpha=0.5)
    axes[2].set_ylabel('Phase (rad)')
    axes[2].set_xlabel('Time (s)')
    axes[2].set_title('Phase Estimation')
    axes[2].legend()
    axes[2].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('/home/zhang/gait_control_ws/test_frequency_change.png', dpi=150)
    print(f"图像已保存: test_frequency_change.png")


def test_noisy_signal():
    """测试3: 带噪声的步态信号"""
    print("\n" + "="*60)
    print("测试3: 带噪声的步态信号")
    print("="*60)
    
    dt = 0.02
    duration = 10.0
    t = np.arange(0, duration, dt)
    n_samples = len(t)
    
    # 生成信号
    gait_frequency = 1.2
    amplitude = 0.3
    noise_level = 0.05  # 噪声幅度
    
    angle_left = generate_gait_signal(t, gait_frequency, amplitude, 0)
    angle_right = generate_gait_signal(t, gait_frequency, amplitude, np.pi)
    
    # 添加高斯噪声
    noise = np.random.normal(0, noise_level, n_samples)
    angle_diff_clean = angle_left - angle_right
    angle_diff_noisy = angle_diff_clean + noise
    
    true_phase = np.mod(2 * np.pi * gait_frequency * t, 2 * np.pi)
    
    # 初始化振荡器
    oscillator = AdaptiveOscillator(initial_frequency=1.0, dt=dt)
    init_idx = 10
    dq_init = (angle_diff_noisy[init_idx] - angle_diff_noisy[init_idx-1]) / dt
    oscillator.initialize(angle_diff_noisy[init_idx], dq_init, gait_frequency)
    
    # 运行振荡器
    estimated_phase = np.zeros(n_samples)
    oscillator_output = np.zeros(n_samples)
    
    for i in range(n_samples):
        if i >= init_idx:
            estimated_phase[i] = oscillator.update(angle_diff_noisy[i], gait_frequency)
            oscillator_output[i] = oscillator.y
    
    # 绘图
    fig, axes = plt.subplots(3, 1, figsize=(12, 10))
    
    # 信号对比
    axes[0].plot(t, angle_diff_clean, 'b-', label='Clean Signal', linewidth=1.5, alpha=0.7)
    axes[0].plot(t, angle_diff_noisy, 'r-', label='Noisy Signal', linewidth=1, alpha=0.5)
    axes[0].plot(t, oscillator_output, 'g--', label='Filtered Output', linewidth=1.5)
    axes[0].set_ylabel('Angle (rad)')
    axes[0].set_title(f'Test 3: Noisy Signal (SNR ≈ {20*np.log10(amplitude/noise_level):.1f} dB)')
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)
    
    # 相位估计
    axes[1].plot(t, true_phase, 'k-', label='True Phase', linewidth=2)
    axes[1].plot(t, estimated_phase, 'b--', label='Estimated Phase', linewidth=1.5)
    axes[1].set_ylabel('Phase (rad)')
    axes[1].set_title('Phase Estimation with Noise')
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)
    
    # 相位误差
    phase_error = np.abs(np.mod(estimated_phase - true_phase + np.pi, 2*np.pi) - np.pi)
    axes[2].plot(t, np.rad2deg(phase_error), 'r-', linewidth=1.5)
    axes[2].set_ylabel('Phase Error (deg)')
    axes[2].set_xlabel('Time (s)')
    axes[2].set_title('Phase Error')
    axes[2].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('/home/zhang/gait_control_ws/test_noisy_signal.png', dpi=150)
    print(f"图像已保存: test_noisy_signal.png")
    
    # 统计信息
    converged_idx = int(2.0 / dt)
    if converged_idx < n_samples:
        mean_error = np.mean(phase_error[converged_idx:])
        print(f"\n收敛后平均相位误差: {np.rad2deg(mean_error):.3f}°")


if __name__ == "__main__":
    print("\n" + "="*60)
    print("2阶自适应振荡器测试程序")
    print("基于 adaptive_oscillators.m 的Python实现")
    print("="*60)
    
    # 运行测试
    test_constant_frequency()
    test_frequency_change()
    test_noisy_signal()
    
    print("\n" + "="*60)
    print("所有测试完成！")
    print("="*60)
    
    plt.show()
