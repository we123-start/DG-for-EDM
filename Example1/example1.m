clear
clc
%% basic value
x_left=0;x_right=1;%区间长度
J=80;%空间剖分
a=1;
Dapp=0.002;%扩散项的系数
beta1=10;%惩罚参数
beta2=1/12;%惩罚参数
h=(x_right-x_left)/J;%剖分长度
% dt=0.01;
%N_t=L/(2*Dapp); %N_t=(L*u)/(2*Dapp);
delta_t=0.001;
dt=delta_t;%时间剖分步长
final_t=0.6;%终止时间
%dt=final_t/N_t;
for j=1:J+1
    x(j)=x_left+(j-1)*h;%区间点
end
p=1;%基函数的阶数
p1=p+1;%积分点个数
[xl,w]=gauleg(p1);%积分点x1以及权重w
[P]=Lpoly(p,xl);%相应积分点的勒让德多项式
[D_P]=D_Lpoly(p,xl);%相应积分点的勒让德多项式的导数
PP=P';
%% the initial value of u
for j=1:J
    x_loc(j,1)=x(j);
    x_loc(j,p1+2)=x(j)+h;
    for k=1:p1
        tl=h/2*(xl(k)+1)+x(j);
        x_loc(j,k+1)=tl;
    end
    for i=1:p1
        if (x_loc(j,i+1)>=0.2)&&(x_loc(j,i+1)<=0.4)
            f(i,j)=sin(pi*(x_loc(j,i+1)-0.2)/0.2);
        else
            f(i,j)=0;
        end
    end
    u(:,1,j)=PP\f(:,j);
end
%% DG method
for s=1:p1
    temp(s)=(-1)^(s-1);
end
temp1=temp';
vxminus=[0;1].*(2/h);
vxplus=[0;1].*(2/h);
%for j=1:J
    for i=0:p
        for m=0:p
            sum_loc1=sum(w.*P(i+1,:).*P(m+1,:));
            M(i+1,m+1)=h/2*sum_loc1;
            sum_loc2=sum(w.*D_P(m+1,:).*D_P(i+1,:));
            D_xx(i+1,m+1)=2/h*sum_loc2;
            sum_loc3=sum(w.*P(m+1,:).*D_P(i+1,:));
            D_x(i+1,m+1)=sum_loc3;
        end
    end
%end
time=0;
kk=1;
while(time<final_t)
    if time+dt>final_t
        dt=final_t-time;
    end
    kk=kk+1;
%     %%%%%% Euler method;
%     for j=1:J
%         for i=0:p%测试函数的个数
%             if j==1
%                 uplus=sum(u(:,kk-1,j+1).*temp1(:));
%                 uminus=sum(u(:,kk-1,j).*ones(p1,1));
%                 penalty=h\beta1*(uplus-uminus);
%                 uxplus=u(2,kk-1,j+1)*2/h+u(3,kk-1,j+1)*2/h*3*(-1)+u(4,kk-1,j+1)*2/h*6;
%                 uxminus=u(2,kk-1,j)*2/h+u(3,kk-1,j)*2/h*3+u(4,kk-1,j)*2/h*6;
%                 uxxplus=u(3,kk-1,j+1)*(2/h)^2*3+u(4,kk-1,j+1)*(2/h)^2*15*(-1);
%                 uxxminus=u(3,kk-1,j)*(2/h)^2*3+u(4,kk-1,j)*(2/h)^2*15;
%                 Ar=Dapp*(penalty+1/2*(uxplus+uxminus)+beta2*h*(uxxplus-uxxminus))-1/2*(uplus+uminus)-Dapp*(1/2*(uplus+uminus)-uminus)*vxminus(i+1);
%                 uplus=sum(u(:,kk-1,j).*temp1(:));
%                 uminus=0;%边界条件
%                 penalty=h\beta1*(uplus-uminus);
%                 uxplus=u(2,kk-1,j)*2/h+u(3,kk-1,j)*2/h*3*(-1)+u(4,kk-1,j)*2/h*6;
%                 uxminus=0;
%                 uxxplus=u(3,kk-1,j)*(2/h)^2*3+u(4,kk-1,j)*(2/h)^2*15*(-1);
%                 uxxminus=0;
%                 Al=(Dapp*(penalty+1/2*(uxplus+uxminus)+beta2*h*(uxxplus-uxxminus))-1/2*(uplus+uminus))*(-1).^i-Dapp*(1/2*(uplus+uminus)-uplus)*vxplus(i+1);
%             elseif j==J
%                 uplus=0;%边界条件
%                 uminus=sum(u(:,kk-1,j).*ones(p1,1));
%                 penalty=h\beta1*(uplus-uminus);
%                 uxplus=0;
%                 uxminus=u(2,kk-1,j)*2/h+u(3,kk-1,j)*2/h*3+u(4,kk-1,j)*2/h*6;
%                 uxxplus=0;
%                 uxxminus=u(3,kk-1,j)*(2/h)^2*3+u(4,kk-1,j)*(2/h)^2*15;
%                 Ar=Dapp*(penalty+1/2*(uxplus+uxminus)+beta2*h*(uxxplus-uxxminus))-1/2*(uplus+uminus)-Dapp*(1/2*(uplus+uminus)-uminus)*vxminus(i+1);
%                 uplus=sum(u(:,kk-1,j).*temp1(:));
%                 uminus=sum(u(:,kk-1,j-1).*ones(p1,1));
%                 penalty=h\beta1*(uplus-uminus);
%                 uxplus=u(2,kk-1,j)*2/h+u(3,kk-1,j)*2/h*3*(-1)+u(4,kk-1,j)*2/h*6;
%                 uxminus=u(2,kk-1,j-1)*2/h+u(3,kk-1,j-1)*2/h*3+u(4,kk-1,j-1)*2/h*6;
%                 uxxplus=u(3,kk-1,j)*(2/h)^2*3+u(4,kk-1,j)*(2/h)^2*15*(-1);
%                 uxxminus=u(3,kk-1,j-1)*(2/h)^2*3+u(4,kk-1,j-1)*(2/h)^2*15;
%                 Al=(Dapp*(penalty+1/2*(uxplus+uxminus)+beta2*h*(uxxplus-uxxminus))-1/2*(uplus+uminus))*(-1).^i-Dapp*(1/2*(uplus+uminus)-uplus)*vxplus(i+1);
%             else
%                 uplus=sum(u(:,kk-1,j+1).*temp1(:));
%                 uminus=sum(u(:,kk-1,j).*ones(p1,1));
%                 penalty=h\beta1*(uplus-uminus);
%                 uxplus=u(2,kk-1,j+1)*2/h+u(3,kk-1,j+1)*2/h*3*(-1)+u(4,kk-1,j+1)*2/h*6;
%                 uxminus=u(2,kk-1,j)*2/h+u(3,kk-1,j)*2/h*3+u(4,kk-1,j)*2/h*6;
%                 uxxplus=u(3,kk-1,j+1)*(2/h)^2*3+u(4,kk-1,j+1)*(2/h)^2*15*(-1);
%                 uxxminus=u(3,kk-1,j)*(2/h)^2*3+u(4,kk-1,j)*(2/h)^2*15;
%                 Ar=Dapp*(penalty+1/2*(uxplus+uxminus)+beta2*h*(uxxplus-uxxminus))-1/2*(uplus+uminus)-Dapp*(1/2*(uplus+uminus)-uminus)*vxminus(i+1);
%                 uplus=sum(u(:,kk-1,j).*temp1(:));
%                 uminus=sum(u(:,kk-1,j-1).*ones(p1,1));
%                 penalty=h\beta1*(uplus-uminus);
%                 uxplus=u(2,kk-1,j)*2/h+u(3,kk-1,j)*2/h*3*(-1)+u(4,kk-1,j)*2/h*6;
%                 uxminus=u(2,kk-1,j-1)*2/h+u(3,kk-1,j-1)*2/h*3+u(4,kk-1,j-1)*2/h*6;
%                 uxxplus=u(3,kk-1,j)*(2/h)^2*3+u(4,kk-1,j)*(2/h)^2*15*(-1);
%                 uxxminus=u(3,kk-1,j-1)*(2/h)^2*3+u(4,kk-1,j-1)*(2/h)^2*15;
%                 Al=(Dapp*(penalty+1/2*(uxplus+uxminus)+beta2*h*(uxxplus-uxxminus))-1/2*(uplus+uminus))*(-1).^i-Dapp*(1/2*(uplus+uminus)-uplus)*vxplus(i+1);
%             end
%             A(i+1)=Ar-Al;
%         end
%         u(1:p1,kk,j)=u(1:p1,kk-1,j)+(2*M)\(D_x*u(1:p1,kk-1,j)-Dapp*D_xx*u(1:p1,kk-1,j)+A')*dt;
%     end
%     %%% the limiter function down
%     for j=1:J
%         if j==1
%             q=p1;
%             while(q>1)
%                 aa=(2*q-3)*u(q,kk-1,j);
%                 b=u(q-1,kk-1,j+1)-u(q-1,kk-1,j);
%                 %                 c=u(q-1,kk,j)-u(q-1,kk,j-1);
%                 %                 c=u(q-1,kk,j)-u(q-1,kk,J);%Periodic boundary conditions
%                 c=u(q-1,kk,j);%left boundary conditions
% %                 if time+dt>0.2  %left boundary conditions
% %                     c=u(q-1,kk-1,j)-0;
% %                 else
% %                     c=u(q-1,kk-1,j)-1;
% %                 end
%                 d=minmod(aa,b,c);
%                 u(q,kk-1,j)=1/(2*q-3)*d;
%                 q=q-1;
%             end
%         elseif j==J
%             q=p1;
%             while(q>1)
%                 aa=(2*q-3)*u(q,kk-1,j);
%                 %                 b=u(q-1,kk,j+1)-u(q-1,kk,j);
%                 %                 b=u(q-1,kk,1)-u(q-1,kk,j);%Periodic boundary conditions
%                 b=-u(q-1,kk-1,j);  %right boundary conditions
%                 c=u(q-1,kk-1,j)-u(q-1,kk-1,j-1);
%                 d=minmod(aa,b,c);
%                 u(q,kk-1,j)=1/(2*q-3)*d;
%                 q=q-1;
%             end
%         else
%             q=p1;
%             while(q>1)
%                 aa=(2*q-3)*u(q,kk-1,j);
%                 b=u(q-1,kk-1,j+1)-u(q-1,kk-1,j);
%                 c=u(q-1,kk-1,j)-u(q-1,kk-1,j-1);
%                 d=minmod(aa,b,c);
%                 u(q,kk-1,j)=1/(2*q-3)*d;
%                 q=q-1;
%             end
%         end
%     end
%     %%% the limiter function up

    %%% 3-Runge-Kutta method;
    % round one
    for j=1:J
        for i=0:p
            if j==1
                uplus=sum(u(:,kk-1,j+1).*temp1(:));
                uminus=sum(u(:,kk-1,j).*ones(p1,1));
                penalty=h\beta1*(uplus-uminus);
                uxplus=u(2,kk-1,j+1)*2/h;
                uxminus=u(2,kk-1,j)*2/h;
                uxxplus=0;
                uxxminus=0;
                Ar=Dapp*(penalty+1/2*(uxplus+uxminus)+beta2*h*(uxxplus-uxxminus))-1/2*(uplus+uminus)-Dapp*(1/2*(uplus+uminus)-uminus)*vxminus(i+1);
                uplus=sum(u(:,kk-1,j).*temp1(:));
                uminus=0;
                penalty=h\beta1*(uplus-uminus);
                uxplus=u(2,kk-1,j)*2/h;
                uxminus=0;
                uxxplus=0;
                uxxminus=0;
                Al=(Dapp*(penalty+1/2*(uxplus+uxminus)+beta2*h*(uxxplus-uxxminus))-1/2*(uplus+uminus))*(-1).^i-Dapp*(1/2*(uplus+uminus)-uplus)*vxplus(i+1);
            elseif j==J
                uplus=sum(u(:,kk-1,J).*temp1(:));
                uminus=sum(u(:,kk-1,j).*ones(p1,1));
                penalty=h\beta1*(uplus-uminus);
                uxplus=u(2,kk-1,J)*2/h;
                uxminus=u(2,kk-1,j)*2/h;
                uxxplus=0;
                uxxminus=0;
                Ar=Dapp*(penalty+1/2*(uxplus+uxminus)+beta2*h*(uxxplus-uxxminus))-1/2*(uplus+uminus)-Dapp*(1/2*(uplus+uminus)-uminus)*vxminus(i+1);
                uplus=sum(u(:,kk-1,j).*temp1(:));
                uminus=sum(u(:,kk-1,j-1).*ones(p1,1));
                penalty=h\beta1*(uplus-uminus);
                uxplus=u(2,kk-1,j)*2/h;
                uxminus=u(2,kk-1,j-1)*2/h;
                uxxplus=0;
                uxxminus=0;
                Al=(Dapp*(penalty+1/2*(uxplus+uxminus)+beta2*h*(uxxplus-uxxminus))-1/2*(uplus+uminus))*(-1).^i-Dapp*(1/2*(uplus+uminus)-uplus)*vxplus(i+1);
            else
                uplus=sum(u(:,kk-1,j+1).*temp1(:));
                uminus=sum(u(:,kk-1,j).*ones(p1,1));
                penalty=h\beta1*(uplus-uminus);
                uxplus=u(2,kk-1,j+1)*2/h;
                uxminus=u(2,kk-1,j)*2/h;
                uxxplus=0;
                uxxminus=0;
                Ar=Dapp*(penalty+1/2*(uxplus+uxminus)+beta2*h*(uxxplus-uxxminus))-1/2*(uplus+uminus)-Dapp*(1/2*(uplus+uminus)-uminus)*vxminus(i+1);
                uplus=sum(u(:,kk-1,j).*temp1(:));
                uminus=sum(u(:,kk-1,j-1).*ones(p1,1));
                penalty=h\beta1*(uplus-uminus);
                uxplus=u(2,kk-1,j)*2/h;
                uxminus=u(2,kk-1,j-1)*2/h;
                uxxplus=0;
                uxxminus=0;
                Al=(Dapp*(penalty+1/2*(uxplus+uxminus)+beta2*h*(uxxplus-uxxminus))-1/2*(uplus+uminus))*(-1).^i-Dapp*(1/2*(uplus+uminus)-uplus)*vxplus(i+1);
            end
            A(i+1)=Ar-Al;
        end
        u1(1:p1,kk-1,j)=u(1:p1,kk-1,j)+((1+a).*M)\(D_x*u(1:p1,kk-1,j)-Dapp*D_xx*u(1:p1,kk-1,j)+A')*dt;
    end
    % round two
    for j=1:J
        for i=0:p
            if j==1
                uplus=sum(u1(:,kk-1,j+1).*temp1(:));
                uminus=sum(u1(:,kk-1,j).*ones(p1,1));
                penalty=h\beta1*(uplus-uminus);
                uxplus=u1(2,kk-1,j+1)*2/h;
                uxminus=u1(2,kk-1,j)*2/h;
                uxxplus=0;
                uxxminus=0;
                Ar=Dapp*(penalty+1/2*(uxplus+uxminus)+beta2*h*(uxxplus-uxxminus))-1/2*(uplus+uminus)-Dapp*(1/2*(uplus+uminus)-uminus)*vxminus(i+1);
                uplus=sum(u1(:,kk-1,j).*temp1(:));
                uminus=0;
                penalty=h\beta1*(uplus-uminus);
                uxplus=u1(2,kk-1,j)*2/h;
                uxminus=0;
                uxxplus=0;
                uxxminus=0;
                Al=(Dapp*(penalty+1/2*(uxplus+uxminus)+beta2*h*(uxxplus-uxxminus))-1/2*(uplus+uminus))*(-1).^i-Dapp*(1/2*(uplus+uminus)-uplus)*vxplus(i+1);
            elseif j==J
                uplus=sum(u(:,kk-1,J).*temp1(:));
                uminus=sum(u(:,kk-1,j).*ones(p1,1));
                penalty=h\beta1*(uplus-uminus);
                uxplus=u(2,kk-1,J)*2/h;
                uxminus=u(2,kk-1,j)*2/h;
                uxxplus=0;
                uxxminus=0;
                Ar=Dapp*(penalty+1/2*(uxplus+uxminus)+beta2*h*(uxxplus-uxxminus))-1/2*(uplus+uminus)-Dapp*(1/2*(uplus+uminus)-uminus)*vxminus(i+1);
                uplus=sum(u1(:,kk-1,j).*temp1(:));
                uminus=sum(u1(:,kk-1,j-1).*ones(p1,1));
                penalty=h\beta1*(uplus-uminus);
                uxplus=u1(2,kk-1,j)*2/h;
                uxminus=u1(2,kk-1,j-1)*2/h;
                uxxplus=0;
                uxxminus=0;
                Al=(Dapp*(penalty+1/2*(uxplus+uxminus)+beta2*h*(uxxplus-uxxminus))-1/2*(uplus+uminus))*(-1).^i-Dapp*(1/2*(uplus+uminus)-uplus)*vxplus(i+1);
            else
                uplus=sum(u1(:,kk-1,j+1).*temp1(:));
                uminus=sum(u1(:,kk-1,j).*ones(p1,1));
                penalty=h\beta1*(uplus-uminus);
                uxplus=u1(2,kk-1,j+1)*2/h;
                uxminus=u1(2,kk-1,j)*2/h;
                uxxplus=0;
                uxxminus=0;
                Ar=Dapp*(penalty+1/2*(uxplus+uxminus)+beta2*h*(uxxplus-uxxminus))-1/2*(uplus+uminus)-Dapp*(1/2*(uplus+uminus)-uminus)*vxminus(i+1);
                uplus=sum(u1(:,kk-1,j).*temp1(:));
                uminus=sum(u1(:,kk-1,j-1).*ones(p1,1));
                penalty=h\beta1*(uplus-uminus);
                uxplus=u1(2,kk-1,j)*2/h;
                uxminus=u1(2,kk-1,j-1)*2/h;
                uxxplus=0;
                uxxminus=0;
                Al=(Dapp*(penalty+1/2*(uxplus+uxminus)+beta2*h*(uxxplus-uxxminus))-1/2*(uplus+uminus))*(-1).^i-Dapp*(1/2*(uplus+uminus)-uplus)*vxplus(i+1);
            end
            A(i+1)=Ar-Al;
        end
        u2(1:p1,kk-1,j)=3/4*u(1:p1,kk-1,j)+1/4*(u1(1:p1,kk-1,j)+((1+a).*M)\(D_x*u1(1:p1,kk-1,j)-Dapp*D_xx*u1(1:p1,kk-1,j)+A')*dt);
    end
    % round three
    for j=1:J
        for i=0:p
            if j==1
                uplus=sum(u2(:,kk-1,j+1).*temp1(:));
                uminus=sum(u2(:,kk-1,j).*ones(p1,1));
                penalty=h\beta1*(uplus-uminus);
                uxplus=u2(2,kk-1,j+1)*2/h;
                uxminus=u2(2,kk-1,j)*2/h;
                uxxplus=0;
                uxxminus=0;
                Ar=Dapp*(penalty+1/2*(uxplus+uxminus)+beta2*h*(uxxplus-uxxminus))-1/2*(uplus+uminus)-Dapp*(1/2*(uplus+uminus)-uminus)*vxminus(i+1);
                uplus=sum(u2(:,kk-1,j).*temp1(:));
                uminus=0;
                penalty=h\beta1*(uplus-uminus);
                uxplus=u2(2,kk-1,j)*2/h;
                uxminus=0;
                uxxplus=0;
                uxxminus=0;
                Al=(Dapp*(penalty+1/2*(uxplus+uxminus)+beta2*h*(uxxplus-uxxminus))-1/2*(uplus+uminus))*(-1).^i-Dapp*(1/2*(uplus+uminus)-uplus)*vxplus(i+1);
            elseif j==J
                uplus=sum(u(:,kk-1,J).*temp1(:));
                uminus=sum(u(:,kk-1,j).*ones(p1,1));
                penalty=h\beta1*(uplus-uminus);
                uxplus=u(2,kk-1,J)*2/h;
                uxminus=u(2,kk-1,j)*2/h;
                uxxplus=0;
                uxxminus=0;
                Ar=Dapp*(penalty+1/2*(uxplus+uxminus)+beta2*h*(uxxplus-uxxminus))-1/2*(uplus+uminus)-Dapp*(1/2*(uplus+uminus)-uminus)*vxminus(i+1);
                uplus=sum(u2(:,kk-1,j).*temp1(:));
                uminus=sum(u2(:,kk-1,j-1).*ones(p1,1));
                penalty=h\beta1*(uplus-uminus);
                uxplus=u2(2,kk-1,j)*2/h;
                uxminus=u2(2,kk-1,j-1)*2/h;
                uxxplus=0;
                uxxminus=0;
                Al=(Dapp*(penalty+1/2*(uxplus+uxminus)+beta2*h*(uxxplus-uxxminus))-1/2*(uplus+uminus))*(-1).^i-Dapp*(1/2*(uplus+uminus)-uplus)*vxplus(i+1);
            else
                uplus=sum(u2(:,kk-1,j+1).*temp1(:));
                uminus=sum(u2(:,kk-1,j).*ones(p1,1));
                penalty=h\beta1*(uplus-uminus);
                uxplus=u2(2,kk-1,j+1)*2/h;
                uxminus=u2(2,kk-1,j)*2/h;
                uxxplus=0;
                uxxminus=0;
                Ar=Dapp*(penalty+1/2*(uxplus+uxminus)+beta2*h*(uxxplus-uxxminus))-1/2*(uplus+uminus)-Dapp*(1/2*(uplus+uminus)-uminus)*vxminus(i+1);
                uplus=sum(u2(:,kk-1,j).*temp1(:));
                uminus=sum(u2(:,kk-1,j-1).*ones(p1,1));
                penalty=h\beta1*(uplus-uminus);
                uxplus=u2(2,kk-1,j)*2/h;
                uxminus=u2(2,kk-1,j-1)*2/h;
                uxxplus=0;
                uxxminus=0;
                Al=(Dapp*(penalty+1/2*(uxplus+uxminus)+beta2*h*(uxxplus-uxxminus))-1/2*(uplus+uminus))*(-1).^i-Dapp*(1/2*(uplus+uminus)-uplus)*vxplus(i+1);
            end
            A(i+1)=Ar-Al;
        end
        u(1:p1,kk,j)=1/3*u(1:p1,kk-1,j)+2/3*(u2(1:p1,kk-1,j)+((1+a).*M)\(D_x*u2(1:p1,kk-1,j)-Dapp*D_xx*u2(1:p1,kk-1,j)+A')*dt);
    end
    time=time+dt;
end
%% return the value of function u
for k=1:kk
    for j=1:J
        unum(k,1+(j-1)*p1:j*p1)=PP*u(:,k,j);
    end
end
%% x-coordinate
for j=1:J
    xx(1+(j-1)*p1:j*p1)=x_loc(j,2:p1+1);
end
%% plot the result
dt=delta_t;
% dt=0.01;
%dt=final_t/N_t;
time=-dt;
kk=0;
while(time<final_t)
    if time+dt>final_t
        dt=final_t-time;
    end
    kk=kk+1;
    time=time+dt;
    cexact(kk,:)=exactfunction(Dapp, xx, time,a);
    plot(xx,unum(kk,:),'r*',xx,cexact(kk,:),'b-')
    axis([x_left,x_right,0,1]);
    legend('numerical solution','exact solution')
    xlabel('x (cm)');
    ylabel('c (g/l)');
    pause(0.00001)
end
vnum=textread('vnum.txt');
vexact=textread('vexact.txt');
plot(xx,unum(end,:),'r*',xx,cexact(end,:),'b-',xx,vnum(end,:),'c*',xx,vexact(end,:),'k-')
legend('numerical solution c1','exact solution c1','numerical solution c2','exact solution c2')
xlabel('x (cm)');
ylabel('c (g/l)');




%  writematrix(unum(end,:), 'vnum.txt', 'Delimiter', ' ');
%  writematrix(cexact(end,:), 'vexact.txt', 'Delimiter', ' ');


% hold on
% cexact(:)=exactfunction(Dapp, xx, final_t);
% plot(xx,cexact(:),'b-')
% err1=max(abs(unum(end,:)-cexact(end,:)));%L_infty
% err2=sqrt((2*J)\sum((unum(end,:)-cexact(end,:)).^2));%L2
% writematrix(err1, '5N1=320.txt', 'Delimiter', ' ');
% writematrix(err2, '5N2=320.txt', 'Delimiter', ' ');

% data1_1=textread('5N1=20.txt');
% data1_2=textread('5N1=40.txt');
% data1_3=textread('5N1=80.txt');
% data1_4=textread('5N1=160.txt');
% data1_5=textread('5N1=320.txt');
% data2_1=textread('5N2=20.txt');
% data2_2=textread('5N2=40.txt');
% data2_3=textread('5N2=80.txt');
% data2_4=textread('5N2=160.txt');
% data2_5=textread('5N2=320.txt');
% cr1_1=log2(data1_1/data1_2);
% cr1_2=log2(data1_2/data1_3);
% cr1_3=log2(data1_3/data1_4);
% cr1_4=log2(data1_4/data1_5);
% cr2_1=log2(data2_1/data2_2);
% cr2_2=log2(data2_2/data2_3);
% cr2_3=log2(data2_3/data2_4);
% cr2_4=log2(data2_4/data2_5);