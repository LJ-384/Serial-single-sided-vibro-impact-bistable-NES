function [eta,eta1,tout,xout]=SSVIT(tfinal)
epsilon1=0.025;  %m2/m1的比值
epsilon2=0.025;  %m3/m1的比值
k1=0.05;       %K2/K1的比值
k2=0.05;       %K3/K1的比值
alpha=0.9;     %L0/L的比值
lambda1=0.001;   %C1/M1*w0
lambda2=0.01; %C2/M1*w0
lambda3=0.01;   %C3/M1*w0
R=sqrt(1-alpha^2);
rc=0.7;     %碰撞损失

%初始条件
v0=2;
%% 运动微分方程(单边)
    function dx=wffc(t,x)   %等价一次微分方程组
        dx=zeros(6,1);
        x1=x(1);v1=x(2);x2=x(3);v2=x(4);x3=x(5);v3=x(6);  %x依次是m1的位移，m1的速度，m2的位移，m2的速度，m3的位移，m3的速度
        Lt=1/sqrt(alpha^2+(x2-x1)^2); 
        dx(1)=v1;  %m1的速度
        dx(2)=(-x1-2*k1*(x1-x2)*(1-Lt)-lambda1*v1-lambda2*(v1-v2)); %m1的加速度
        dx(3)=v2;  %m2的速度
        dx(4)=(1/epsilon1)*(-k2*(x2-x3)-2*k1*(x2-x1)*(1-Lt)-lambda2*(v2-v1)-lambda3*(v2-v3));  %m2的加速度
        dx(5)=v3; %m3的速度
        dx(6)=(1/epsilon2)*(-k2*(x3-x2)-lambda3*(v3-v2)); %m3的加速度
        if snap==1    %m3和m2在正方向黏住
            if dx(4)>dx(2)
                dx(2)=-(k1*(x2-x1)+c2*(v2-v1))/(m1+m2);
                dx(4)=dx(2);
            else
                snap=0;
            end
        end
    end
%% 碰撞检测函数，求解运动微分方程时调用该函数检测是否发生碰撞
    function [value,isterminal,direction]=events(t,x)
        value=x(1)-x(3)-R;     % detect gap = 0
        isterminal = 1;   % stop the integration
        direction = 1;   % positive direction
    end
%% 求解
    x0=[0 v0 -R 0 -R 0];                      % 0初始条件,M1和M2的初始速度为v0，m3初始位置在zc处，此时弹簧无压缩
    tstart = 0;     %求解起始时间
    options = odeset('Events',@events); %求解终断判断设置，调用上述的碰撞检测函数，当间隙为zc时停止计算
    tout = tstart;  %用于存储最终输出的时间列
    xout = x0;      %用于存储最终输出的各个自由度位移与速度。
    teout = [];     %用于存储发生碰撞的时间
    xeout = [];     %用于存储发生碰撞时各个自由度的位移和速度
    ieout = [];     %用于存储发生碰撞的event，这里只设置了一个event，所以这个变量没什么用。
    snap=0;         %默认m3和m2没有黏住
    while tstart<tfinal %循环计算，直到计算到所需的时间
       [t,x,te,xe,ie]=ode45(@wffc,[tstart tfinal],x0,options);      %求解振动微分方程，发生碰撞后中断
       nt = length(t);
       tout = [tout; t(2:nt)];              %将求解得到的时间列存入tou
       xout = [xout; x(2:nt,:)];            %将求解得到的各个自由度位移速度存入xout
       % Set the new initial conditions, with attenuation.
       if ~isempty(ie) % 如果ie不为空，表示发生了碰撞而计算终断
           teout = [teout; te];             % 将碰撞时的时间存入teout
           xeout = [xeout; xe];             % 将碰撞时的各个自由度位移速度存入xeout
           ieout = [ieout; ie];
           x0=xe;                           %发生碰撞后，根据公式计算碰撞后的速度，然后将碰撞后的状态作为初始状态，再次求解。此处x0为下次求解的初始状态。
           v2m=xe(4);                                   %碰撞前m2的速度；
           v1m=xe(2);                                   %碰撞前m1的速度；
           v1p=(v1m*(1-epsilon1*rc)+v2m*epsilon1*(1+rc))/(1+epsilon1);    %碰撞后m1的速度；
           v2p=(v1m*(1+rc)+v2m*(epsilon1-rc))/(1+epsilon1);         %碰撞后m2的速度；
           x0(4)=v2p;                     %将初始条件中m2和m3的速度设置为碰撞后的速度
           x0(2)=v1p;                     %将初始条件中m2和m3的速度设置为碰撞后的速度
           %options = odeset(options,'InitialStep',t(nt)-t(nt-1)); %将上一段的计算步长作为下一次计算的初始步长
        if x0(3)<x0(1)&&x0(4)-x0(2)<=1e-5           %负向碰撞：如果发生碰撞后m2和m3的速度差非常小，可以认为两者速度一致，由于碰撞时位移也一致，视为粘在一起了，修改snap的值。
            snap=2;
        end
       end
       tstart = t(nt);    %改变下一次计算的初始时间为上一段计算的结束时间。
    end
%% 能量耗散
    E0=0.5*v0^2;
    eta_impact=0;
    if ~isempty(teout)  %如果发生了碰撞
        v2m=xeout(:,4); %提取每次碰撞前的速度
        v1m=xeout(:,2);
        v1p=(v1m*(1-epsilon1*rc)+v2m*epsilon1*(1+rc))/(1+epsilon1); %计算碰撞后的速度
        v2p=(v1m*(1+rc)+v2m*(epsilon1-rc))/(1+epsilon1);
        DKEm1=0.5*(v1m.^2-v1p.^2)/E0;  %m1碰撞前后动能变化
        DKEm2=0.5*epsilon1*(v2m.^2-v2p.^2)/E0;  %m2碰撞前后动能变化
        eta_impact=cumsum(DKEm1+DKEm2);   %两者累加，得到因为碰撞耗散掉的能量
    end
    eta_impact=interp1([teout;0],[eta_impact;0],tout,'previous','extrap'); 
    eta=eta_impact+lambda1*cumtrapz(tout,(xout(:,2)).^2)/E0+lambda2*cumtrapz(tout,(xout(:,4)-xout(:,2)).^2)/E0+lambda3*cumtrapz(tout,(xout(:,6)-xout(:,4)).^2)/E0;
    eta1=0.5*(xout(:,1).^2)/E0+0.5*(xout(:,2).^2)/E0;
end