function [PF,m]=lmd(x)

% 局域均值分析
c = x;
N = length(x);
 
A = ones(1,N);
PF = [];
aii = 2*A;
 
while(1)
 
  si = c;
  a = 1;
  
   while(1)
    h = si;
    
      maxVec = [];
      minVec = [];
      
   % 寻找极大值点和极小值点
      for i = 2: N - 1
         if h (i - 1) < h (i) & h (i) > h (i + 1)
            maxVec = [maxVec i]; 		
         end
         if h (i - 1) > h (i) & h (i) < h (i + 1)
            minVec = [minVec i]; 		
         end         
      end
      
   % 检查是否有残余
      if (length (maxVec) + length (minVec)) < 2
         break;
      end
           
  % 原始信号中的两边两个点的判断 
      lenmax=length(maxVec);
      lenmin=length(minVec);
      %先是左边这个点
      if h(1)>0
          if(maxVec(1)<minVec(1))
              yleft_max=h(maxVec(1));
              yleft_min=-h(1);
          else
              yleft_max=h(1);
              yleft_min=h(minVec(1));
          end
      else
          if (maxVec(1)<minVec(1))
              yleft_max=h(maxVec(1));
              yleft_min=h(1);
          else
              yleft_max=-h(1);
              yleft_min=h(minVec(1));
          end
      end
      %然后判断右边这个点
      if h(N)>0
          if(maxVec(lenmax)<minVec(lenmin))
             yright_max=h(N);
             yright_min=h(minVec(lenmin));
          else
              yright_max=h(maxVec(lenmax));
              yright_min=-h(N);
          end
      else
          if(maxVec(lenmax)<minVec(lenmin))
              yright_max=-h(N);
              yright_min=h(minVec(lenmin));
          else
              yright_max=h(maxVec(lenmax));
              yright_min=h(N);
          end
      end
      %使用三次样条插值方法，对极大值向量和极小值向量进行插值
      %spline interpolate
      maxEnv=spline([1 maxVec N],[yleft_max h(maxVec) yright_max],1:N);
      minEnv=spline([1 minVec N],[yleft_min h(minVec) yright_min],1:N);
      
    mm = (maxEnv + minEnv)/2;%得到局部均值函数
    aa = abs(maxEnv - minEnv)/2;%得到包络函数
    
    mmm = mm;
    aaa = aa;
 
    preh = h;
    h = h-mmm;%从原始信号中分离处局部均值函数
    si = h./aaa;%对分离出的信号进行解调
    a = a.*aaa;    
    
aii = aaa;
 
    B = length(aii);
    C = ones(1,B);
    bb = norm(aii-C);%返回aii-C的最大奇异值，aii就是那个包络函数
    if(bb < 1000)%如果bb<1000，就得到了纯调频函数
        break;
    end     
    
   end %分解1个Pf分量在这结束
   
  pf = a.*si;%包络函数和纯调频函数相乘，得到PF分量
  
  PF = [PF; pf];
  
  bbb = length (maxVec) + length (minVec);
 % 简单的一个结束分解的条件
      if (length (maxVec) + length (minVec)) < 20
         break;
      end
           
  c = c-pf;
 
end
m=x-PF(1,:)-PF(2,:)-PF(3,:);%如果分解出2个，从原始信号中减去，得到残余分量
end