(() => {
 const source=document.getElementById('paired-terms')||document.getElementById('article-terms'),root=source?.id==='paired-terms'?document.querySelector('.paired-article'):document.getElementById('article-content');if(!root||!source)return;
 const terms=JSON.parse(source.textContent),svgLabels=new Map();
 Object.assign(terms,{
  '时间 t / ms':'Time t / ms','积分电位 h':'Integrated potential h','新输入 b':'New input b','LSEε(a, b) / 积分值':'LSEε(a, b) / integrated value','固定 a = 0.7，改变 b':'Fix a = 0.7; vary b','ds/dh / 真实局部导数':'ds/dh / true local derivative','峰值 = 1 / (4ε)':'Peak = 1 / (4ε)','平滑前向 sigmoid':'Smooth forward sigmoid','硬阈值 H':'Hard threshold H','最大间隙 ε ln 2':'Maximum gap ε ln 2','阈值处的峰值':'Peak at threshold','回溯步数 k':'Backward steps k','j < 1：越早的状态，影响越弱。':'j < 1: earlier states have less influence.','j > 1：连乘会逐步放大扰动。':'j > 1: repeated products amplify perturbations.','j = 1：这个标量例子的梯度保持不变。':'j = 1: the scalar gradient stays constant.','减弱 LTD':'Depression LTD','增强 LTP':'Potentiation LTP','前 PRE':'Pre','后 POST':'Post','突触前先放：增强这条突触。':'Pre fires first: potentiate the synapse.','突触后先放：减弱这条突触。':'Post fires first: depress the synapse.','同一时刻：本实验约定不更新。':'Simultaneous: no update by convention.','目标时刻':'Target time','切触：速度 = 0':'Tangency: speed = 0','峰值低于阈值':'Peak below threshold','临界点':'Critical point','无事件':'No event','未定义':'Undefined','w = 4：峰值处 du/dt = 0。':'w = 4: du/dt = 0 at the peak.','w < 4：电压达不到阈值。':'w < 4: voltage does not reach threshold.','此处不定义首次脉冲时刻损失。':'First-spike time loss is undefined here.','阈值被切触而非横截穿越：不能使用事件时间梯度公式。':'Tangency, not a transversal crossing: the event-time derivative does not apply.','该脉冲不存在：不把缺失事件误当作零梯度。':'The spike does not exist: a missing event is not a zero gradient.','已经达到目标时刻，当前梯度为零。':'The target time is reached; the gradient is zero.','梯度下降将增大权重，使放电提前。':'Gradient descent increases the weight, advancing the spike.','J ≤ θ：不能在有限时间内穿越阈值。':'J ≤ θ: no threshold crossing in finite time.','已复制 ✓':'Copied ✓','请选择代码复制':'Select the code to copy','展开全部推导 ＋':'Expand all derivations +','收起全部推导 −':'Collapse all derivations −',
 });
 function english(zh){
  if(terms[zh])return terms[zh];
  if(/^200 ms 内：/.test(zh))return zh.replace('200 ms 内：','Within 200 ms: ').replace(' 个脉冲',' spikes');
  if(/^积分到阈值：/.test(zh))return zh.replace('积分到阈值：','Time to threshold: ');
  if(/^相邻脉冲间隔：/.test(zh))return zh.replace('相邻脉冲间隔：','Interspike interval: ');
  if(zh.includes('ms / 权重单位'))return zh.replace('ms / 权重单位','ms / weight unit');
  if(/^解析 /.test(zh))return zh.replace('解析 ','Analytic ');
  if(zh.includes('中心差分'))return zh.replace('中心差分','Central difference').replace('不适用','Not applicable');
  if(zh.includes('绝对误差'))return zh.replace('绝对误差','Absolute error');
  return null;
 }
 function formatSvg(element,en,zh){const language=document.documentElement.dataset.language;const text=language==='zh'?zh:language==='en'?en:`${en} / ${zh}`;if(element.textContent!==text)element.textContent=text;}
 function translate(node){
  if(node.nodeType===Node.TEXT_NODE){
   const parent=node.parentElement;if(!parent||parent.closest('[data-lang],script,style,code,pre,textarea,option,noscript,.runtime-pair')||!/[\u3400-\u9fff]/.test(node.data))return;
   const zh=node.data.trim(),en=english(zh);if(!en)return;
   if(parent.namespaceURI==='http://www.w3.org/2000/svg'){
    if(parent.tagName==='text'){svgLabels.set(parent,{en,zh});formatSvg(parent,en,zh);}return;
   }
   const wrapper=document.createElement('span');wrapper.className='i18n runtime-pair';for(const [lang,text]of [['en',en],['zh',zh]]){const span=document.createElement('span');span.dataset.lang=lang;span.lang=lang==='zh'?'zh-CN':'en';span.textContent=text;wrapper.append(span);}node.replaceWith(wrapper);return;
  }
  if(node.nodeType===Node.ELEMENT_NODE&&node.matches('[data-lang],script,style,code,pre,textarea,option,noscript,.runtime-pair'))return;
  [...node.childNodes].forEach(translate);
 }
 translate(root);
 const observer=new MutationObserver(records=>{observer.disconnect();for(const r of records){if(r.type==='characterData')translate(r.target);else for(const n of r.addedNodes)translate(n);}for(const node of svgLabels.keys())if(!node.isConnected)svgLabels.delete(node);observer.observe(root,{subtree:true,childList:true,characterData:true});});
 observer.observe(root,{subtree:true,childList:true,characterData:true});
 window.addEventListener('languagechange',()=>{for(const [node,label]of svgLabels){if(!node.isConnected){svgLabels.delete(node);continue;}formatSvg(node,label.en,label.zh);}});
})();
