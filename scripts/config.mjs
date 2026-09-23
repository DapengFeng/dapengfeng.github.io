export const site = {
  name: 'FENG / 知识实验室', author: 'Dapeng Feng', authorZh: '冯大鹏', url: 'https://dapengfeng.github.io',
  description: '用图形理解数学，用实验认识物理，用数据探索代码的性能。冯大鹏的个人知识实验室。',
};
export const categories = [
  { id: 'math', name: '数学与算法', en: 'MATHEMATICS', symbol: '∑', description: '从抽象公式，到直觉与图形', color: '#bcf46e' },
  { id: 'physics', name: '物理与模型', en: 'PHYSICS', symbol: '∿', description: '用模型，解释世界如何运转', color: '#88bfff' },
  { id: 'systems', name: '系统与性能', en: 'SYSTEMS', symbol: '⌘', description: '走近底层，让每个周期有价值', color: '#bca4ff' },
  { id: 'benchmark', name: '基准测试', en: 'EXPERIMENTS', symbol: '▥', description: '少一点猜测，多一点测量', color: '#f3ba83' },
];
export const escape = value => String(value ?? '').replace(/[&<>"']/g, c => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
