export const site = {
  name: 'FENG / 知识实验室', author: 'Dapeng Feng', authorZh: '冯大鹏', url: 'https://dapengfeng.github.io',
  description: '从具体的问题与观察出发，把复杂的原理讲清楚，把值得留意的细节记录下来。冯大鹏的个人知识实验室。',
};
export const categories = [
  { id: 'math', name: '数学与算法', en: 'MATHEMATICS', symbol: '∑', description: '从抽象公式，到直觉与图形', color: '#bcf46e' },
  { id: 'physics', name: '物理与模型', en: 'PHYSICS', symbol: '∿', description: '用模型，解释世界如何运转', color: '#88bfff' },
  { id: 'systems', name: '系统与性能', en: 'SYSTEMS', symbol: '⌘', description: '走近底层，让每个周期有价值', color: '#bca4ff' },
  { id: 'benchmark', name: '基准测试', en: 'EXPERIMENTS', symbol: '▥', description: '少一点猜测，多一点测量', color: '#f3ba83' },
  { id: 'biology', name: '生命与神经科学', en: 'LIFE & NEUROSCIENCE', symbol: '◎', description: '从细胞结构到神经计算', color: '#f2a7bc' },
  { id: 'travel', name: '旅行与见闻', en: 'TRAVEL', symbol: '⌁', description: '沿途的风景，与日常之外的观察', color: '#79d5c6' },
];
export const escape = value => String(value ?? '').replace(/[&<>"']/g, c => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
