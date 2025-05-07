function f = objectivefcn11(N,quantizeda,m_shengcheng,sensors,x)

for i=1:1:N

    AA(i)=quantizeda(1,i)-m_shengcheng(i,1)-sqrt((sensors(1,i)-x(1)).^2+(sensors(2,i)-x(2)).^2);

end

f = AA*AA';
