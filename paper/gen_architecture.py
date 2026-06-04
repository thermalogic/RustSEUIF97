import matplotlib.pyplot as plt
import matplotlib.patches as patches

fig, ax = plt.subplots(figsize=(8, 6))
ax.set_xlim(0, 10)
ax.set_ylim(0, 8)
ax.axis('off')

# Colors
layer_colors = ['#E8F4FD', '#E8F8E8', '#FFF3E0']
border_colors = ['#2196F3', '#4CAF50', '#FF9800']

# Layer 1: Public API Layer
rect1 = patches.FancyBboxPatch((1.5, 5.8), 7, 1.2, boxstyle="round,pad=0.1", 
                                facecolor=layer_colors[0], edgecolor=border_colors[0], linewidth=2)
ax.add_patch(rect1)
ax.text(5, 6.7, 'Public API Layer', ha='center', va='center', fontsize=12, fontweight='bold')
ax.text(5, 6.1, 'pt(p,T,OH) | ph(p,h,OV) | ps(p,s,OD) | pv(p,v,OT) | ...', 
        ha='center', va='center', fontsize=10, family='monospace')

# Arrow 1 (from Layer 1 bottom to Layer 2 top)
ax.annotate('', xy=(5, 5.0), xytext=(5, 5.8),
            arrowprops=dict(arrowstyle='->', color='#333', lw=2))

# Layer 2: Region Wrapper Layer
rect2 = patches.FancyBboxPatch((1.5, 3.8), 7, 1.2, boxstyle="round,pad=0.1",
                                facecolor=layer_colors[1], edgecolor=border_colors[1], linewidth=2)
ax.add_patch(rect2)
ax.text(5, 4.7, 'Region Wrapper Layer', ha='center', va='center', fontsize=12, fontweight='bold')
ax.text(5, 4.1, 'pT_reg1(π,τ) | pT_reg2(π,τ) | pT_reg3(ρ,T) | pT_reg5(π,τ)',
        ha='center', va='center', fontsize=10, family='monospace')

# Arrow 2 (from Layer 2 bottom to Layer 3 top)
ax.annotate('', xy=(5, 3.0), xytext=(5, 3.8),
            arrowprops=dict(arrowstyle='->', color='#333', lw=2))

# Layer 3: Optimized Kernel Layer
rect3 = patches.FancyBboxPatch((1.5, 1.8), 7, 1.2, boxstyle="round,pad=0.1",
                                facecolor=layer_colors[2], edgecolor=border_colors[2], linewidth=2)
ax.add_patch(rect3)
ax.text(5, 2.7, 'Optimized Kernel Layer', ha='center', va='center', fontsize=12, fontweight='bold')
ax.text(5, 2.1, 'poly_powi_steps | poly_j_powi_steps | polys_i_j_powi_steps',
        ha='center', va='center', fontsize=10, family='monospace')

# Labels on the left
ax.text(0.8, 6.4, 'User Interface', ha='center', va='center', fontsize=9, style='italic', color='#555')
ax.text(0.8, 4.4, 'Region Dispatch', ha='center', va='center', fontsize=9, style='italic', color='#555')
ax.text(0.8, 2.4, 'Core Algorithms', ha='center', va='center', fontsize=9, style='italic', color='#555')

plt.tight_layout()
plt.savefig('d:/IF97/IF97-Rust/RustSEUIF97/paper/architecture.png', dpi=300, bbox_inches='tight', facecolor='white')
print("Figure saved successfully")
