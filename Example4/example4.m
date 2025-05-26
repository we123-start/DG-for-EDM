clear
%% basic value
x_left=0;x_right=5;

Nt=5000;
J=200;
L=5;
%Dapp=0.002;
%Dapp=(L*u)/(2*Nt);
Dapp=(L)/(2*Nt);%u=1
beta1=10;
beta2=1/12;
h=(x_right-x_left)/J;
initial_t=0;
final_t=25;
deltat=final_t/Nt;
dt=deltat;
%dt=final_t/Nt;
umin=0.0317;    %picture 4
vmin=0.3016;

f0=0.2051;    %picture 4
g0=0.0598;

for j=1:J+1
    x(j)=x_left+(j-1)*h;
end
p=3;%基函数的阶数
p1=p+1;%积分点个数
[xl,w]=gauleg(p1);%积分点x1以及权重w
[P]=Lpoly(p,xl);%相应积分点的勒让德多项式
[D_P]=D_Lpoly(p,xl);%相应积分点的勒让德多项式的导数
PP=P';
for s=1:p1
    temp(s)=(-1)^(s-1);
end
temp1=temp';
wxminus=[0;1;3;6].*(2/h);
wxplus=[0;1;-3;6].*(2/h);
%% the initial value of u
for j=1:J
    x_loc(j,1)=x(j);
    x_loc(j,p1+2)=x(j)+h;
    for k=1:p1
        tl=h/2*(xl(k)+1)+x(j);
        x_loc(j,k+1)=tl;
    end
    for i=1:p1
        f(i,j)=f0;
        g(i,j)=g0;
    end
    u(:,1,j)=PP\f(:,j);
    v(:,1,j)=PP\g(:,j);
end
%% DG method with forword Euler
time=0;
kk=1;
while(time<final_t)
    if time+dt>final_t
        dt=final_t-time;
    end
    kk=kk+1;
%     %%%    forward Euler
%     for j=1:J
%         for i=0:p
%             for m=0:p
%                 sum_loc1=sum(w.*P(i+1,:).*P(m+1,:));
%                 M(i+1,m+1)=h/2*sum_loc1;
%                 sum_loc2=sum(w.*D_P(m+1,:).*D_P(i+1,:));
%                 D_xx(i+1,m+1)=2/h*sum_loc2;
%                 sum_loc3=sum(w.*P(m+1,:).*D_P(i+1,:));
%                 D_x(i+1,m+1)=sum_loc3;
%                 sum_loc1_1=0;
%                 sum_loc1_2=0;
%                 sum_loc1_3=0;
%                 sum_loc1_4=0;
%                 sum_loc1_5=0;
%                 sum_loc1_6=0;
%                 for k=0:p
%                     sum_loc1_1=sum_loc1_1+w(k+1).*(2-sin(sum(u(1:p1,kk-1,j).*P(1:p+1,k+1)))).*P(i+1,k+1).*P(m+1,k+1);
%                     sum_loc1_2=sum_loc1_2+w(k+1).*0.1.*P(i+1,k+1).*P(m+1,k+1);
%                     %             sum_loc1_3=sum_loc1_3+w(k+1).*(-1.5*sum(u(1:p1,kk-1,j).*P(1:p+1,k+1))/(1-sum(u(1:p1,kk-1,j).*P(1:p+1,k+1))+sum(v(1:p1,kk-1,j).*P(1:p+1,k+1)))^2.)*P(i+1,k+1);
%                     sum_loc1_4=sum_loc1_4+w(k+1).*0.1.*P(i+1,k+1).*P(m+1,k+1);
%                     sum_loc1_5=sum_loc1_5+w(k+1).*(2+cos(sum(v(1:p1,kk-1,j).*P(1:p+1,k+1))))*P(i+1,k+1).*P(m+1,k+1);
%                     %             sum_loc1_6=sum_loc1_6+w(k+1).*(-3*sum(v(1:p1,kk-1,j).*P(1:p+1,k+1))/(1-sum(u(1:p1,kk-1,j).*P(1:p+1,k+1))+sum(v(1:p1,kk-1,j).*P(1:p+1,k+1)))^2.)*P(i+1,k+1);
%                 end
%                 E1(i+1,m+1)=h/2*sum_loc1_1;
%                 F1(i+1,m+1)=h/2*sum_loc1_2;
%                 %                 G1(i+1)=h/2*sum_loc1_3;
%                 E2(i+1,m+1)=h/2*sum_loc1_4;
%                 F2(i+1,m+1)=h/2*sum_loc1_5;
%                 %                 G2(i+1)=h/2*sum_loc1_6;
%             end
%             if j==1
%                 uplus=sum(u(:,kk-1,j+1).*temp1(:));
%                 uminus=sum(u(:,kk-1,j).*ones(p1,1));
%                 penalty=h\beta1*(uplus-uminus);
%                 uxplus=u(2,kk-1,j+1)*2/h+u(3,kk-1,j+1)*2/h*3*(-1)+u(4,kk-1,j+1)*2/h*6;
%                 uxminus=u(2,kk-1,j)*2/h+u(3,kk-1,j)*2/h*3+u(4,kk-1,j)*2/h*6;
%                 uxxplus=u(3,kk-1,j+1)*(2/h)^2*3+u(4,kk-1,j+1)*(2/h)^2*15*(-1);
%                 uxxminus=u(3,kk-1,j)*(2/h)^2*3+u(4,kk-1,j)*(2/h)^2*15;
%                 A1r=Dapp*(penalty+1/2*(uxplus+uxminus)+beta2*h*(uxxplus-uxxminus))-1/2*(uplus+uminus)-Dapp*(1/2*(uplus+uminus)-uminus)*wxminus(i+1);
%                 uplus=sum(u(:,kk-1,j).*temp1(:));
%                 uminus=umin;%left boundary conditions
%                 penalty=h\beta1*(uplus-uminus);
%                 uxplus=u(2,kk-1,j)*2/h+u(3,kk-1,j)*2/h*3*(-1)+u(4,kk-1,j)*2/h*6;
%                 uxminus=0;
%                 uxxplus=u(3,kk-1,j)*(2/h)^2*3+u(4,kk-1,j)*(2/h)^2*15*(-1);
%                 uxxminus=0;
%                 A1l=(Dapp*(penalty+1/2*(uxplus+uxminus)+beta2*h*(uxxplus-uxxminus))-1/2*(uplus+uminus))*(-1).^i-Dapp*(1/2*(uplus+uminus)-uplus)*wxplus(i+1);
%                 vplus=sum(v(:,kk-1,j+1).*temp1(:));
%                 vminus=sum(v(:,kk-1,j).*ones(p1,1));
%                 penalty=h\beta1*(vplus-vminus);
%                 vxplus=v(2,kk-1,j+1)*2/h+v(3,kk-1,j+1)*2/h*3*(-1)+v(4,kk-1,j+1)*2/h*6;
%                 vxminus=v(2,kk-1,j)*2/h+v(3,kk-1,j)*2/h*3+v(4,kk-1,j)*2/h*6;
%                 vxxplus=v(3,kk-1,j+1)*(2/h)^2*3+v(4,kk-1,j+1)*(2/h)^2*15*(-1);
%                 vxxminus=v(3,kk-1,j)*(2/h)^2*3+v(4,kk-1,j)*(2/h)^2*15;
%                 A2r=Dapp*(penalty+1/2*(vxplus+vxminus)+beta2*h*(vxxplus-vxxminus))-1/2*(vplus+vminus)-Dapp*(1/2*(vplus+vminus)-vminus)*wxminus(i+1);
%                 vplus=sum(v(:,kk-1,j).*temp1(:));
%                 vminus=vmin;%left boundary conditions
%                 penalty=h\beta1*(vplus-vminus);
%                 vxplus=v(2,kk-1,j)*2/h+v(3,kk-1,j)*2/h*3*(-1)+v(4,kk-1,j)*2/h*6;
%                 vxminus=0;
%                 vxxplus=v(3,kk-1,j)*(2/h)^2*3+v(4,kk-1,j)*(2/h)^2*15*(-1);
%                 vxxminus=0;
%                 A2l=(Dapp*(penalty+1/2*(vxplus+vxminus)+beta2*h*(vxxplus-vxxminus))-1/2*(vplus+vminus))*(-1).^i-Dapp*(1/2*(vplus+vminus)-vplus)*wxplus(i+1);
%             elseif j==J
%                 uplus=sum(u(:,kk-1,j).*ones(p1,1));%right boundary conditions
%                 uminus=sum(u(:,kk-1,j).*ones(p1,1));
%                 penalty=h\beta1*(uplus-uminus);
%                 uxplus=u(2,kk-1,j)*2/h+u(3,kk-1,j)*2/h*3*(-1)+u(4,kk-1,j)*2/h*6;
%                 uxminus=u(2,kk-1,j)*2/h+u(3,kk-1,j)*2/h*3+u(4,kk-1,j)*2/h*6;
%                 uxxplus=u(3,kk-1,j)*(2/h)^2*3+u(4,kk-1,j)*(2/h)^2*15*(-1);
%                 uxxminus=u(3,kk-1,j)*(2/h)^2*3+u(4,kk-1,j)*(2/h)^2*15;
%                 A1r=Dapp*(penalty+1/2*(uxplus+uxminus)+beta2*h*(uxxplus-uxxminus))-1/2*(uplus+uminus)-Dapp*(1/2*(uplus+uminus)-uminus)*wxminus(i+1);
%                 uplus=sum(u(:,kk-1,j).*temp1(:));
%                 uminus=sum(u(:,kk-1,j-1).*ones(p1,1));
%                 penalty=h\beta1*(uplus-uminus);
%                 uxplus=u(2,kk-1,j)*2/h+u(3,kk-1,j)*2/h*3*(-1)+u(4,kk-1,j)*2/h*6;
%                 uxminus=u(2,kk-1,j-1)*2/h+u(3,kk-1,j-1)*2/h*3+u(4,kk-1,j-1)*2/h*6;
%                 uxxplus=u(3,kk-1,j)*(2/h)^2*3+u(4,kk-1,j)*(2/h)^2*15*(-1);
%                 uxxminus=u(3,kk-1,j-1)*(2/h)^2*3+u(4,kk-1,j-1)*(2/h)^2*15;
%                 A1l=(Dapp*(penalty+1/2*(uxplus+uxminus)+beta2*h*(uxxplus-uxxminus))-1/2*(uplus+uminus))*(-1).^i-Dapp*(1/2*(uplus+uminus)-uplus)*wxplus(i+1);
%                 vplus=sum(v(:,kk-1,j).*ones(p1,1));%right boundary conditions
%                 vminus=sum(v(:,kk-1,j).*ones(p1,1));
%                 penalty=h\beta1*(vplus-vminus);
%                 vxplus=v(2,kk-1,j)*2/h+v(3,kk-1,j)*2/h*3*(-1)+v(4,kk-1,j)*2/h*6;
%                 vxminus=v(2,kk-1,j)*2/h+v(3,kk-1,j)*2/h*3+v(4,kk-1,j)*2/h*6;
%                 vxxplus=v(3,kk-1,j)*(2/h)^2*3+v(4,kk-1,j)*(2/h)^2*15*(-1);
%                 vxxminus=v(3,kk-1,j)*(2/h)^2*3+v(4,kk-1,j)*(2/h)^2*15;
%                 A2r=Dapp*(penalty+1/2*(vxplus+vxminus)+beta2*h*(vxxplus-vxxminus))-1/2*(vplus+vminus)-Dapp*(1/2*(vplus+vminus)-vminus)*wxminus(i+1);
%                 vplus=sum(v(:,kk-1,j).*temp1(:));
%                 vminus=sum(v(:,kk-1,j-1).*ones(p1,1));
%                 penalty=h\beta1*(vplus-vminus);
%                 vxplus=v(2,kk-1,j)*2/h+v(3,kk-1,j)*2/h*3*(-1)+v(4,kk-1,j)*2/h*6;
%                 vxminus=v(2,kk-1,j-1)*2/h+v(3,kk-1,j-1)*2/h*3+v(4,kk-1,j-1)*2/h*6;
%                 vxxplus=v(3,kk-1,j)*(2/h)^2*3+v(4,kk-1,j)*(2/h)^2*15*(-1);
%                 vxxminus=v(3,kk-1,j-1)*(2/h)^2*3+v(4,kk-1,j-1)*(2/h)^2*15;
%                 A2l=(Dapp*(penalty+1/2*(vxplus+vxminus)+beta2*h*(vxxplus-vxxminus))-1/2*(vplus+vminus))*(-1).^i-Dapp*(1/2*(vplus+vminus)-vplus)*wxplus(i+1);
%             else
%                 uplus=sum(u(:,kk-1,j+1).*temp1(:));
%                 uminus=sum(u(:,kk-1,j).*ones(p1,1));
%                 penalty=h\beta1*(uplus-uminus);
%                 uxplus=u(2,kk-1,j+1)*2/h+u(3,kk-1,j+1)*2/h*3*(-1)+u(4,kk-1,j+1)*2/h*6;
%                 uxminus=u(2,kk-1,j)*2/h+u(3,kk-1,j)*2/h*3+u(4,kk-1,j)*2/h*6;
%                 uxxplus=u(3,kk-1,j+1)*(2/h)^2*3+u(4,kk-1,j+1)*(2/h)^2*15*(-1);
%                 uxxminus=u(3,kk-1,j)*(2/h)^2*3+u(4,kk-1,j)*(2/h)^2*15;
%                 A1r=Dapp*(penalty+1/2*(uxplus+uxminus)+beta2*h*(uxxplus-uxxminus))-1/2*(uplus+uminus)-Dapp*(1/2*(uplus+uminus)-uminus)*wxminus(i+1);
%                 uplus=sum(u(:,kk-1,j).*temp1(:));
%                 uminus=sum(u(:,kk-1,j-1).*ones(p1,1));
%                 penalty=h\beta1*(uplus-uminus);
%                 uxplus=u(2,kk-1,j)*2/h+u(3,kk-1,j)*2/h*3*(-1)+u(4,kk-1,j)*2/h*6;
%                 uxminus=u(2,kk-1,j-1)*2/h+u(3,kk-1,j-1)*2/h*3+u(4,kk-1,j-1)*2/h*6;
%                 uxxplus=u(3,kk-1,j)*(2/h)^2*3+u(4,kk-1,j)*(2/h)^2*15*(-1);
%                 uxxminus=u(3,kk-1,j-1)*(2/h)^2*3+u(4,kk-1,j-1)*(2/h)^2*15;
%                 A1l=(Dapp*(penalty+1/2*(uxplus+uxminus)+beta2*h*(uxxplus-uxxminus))-1/2*(uplus+uminus))*(-1).^i-Dapp*(1/2*(uplus+uminus)-uplus)*wxplus(i+1);
%                 vplus=sum(v(:,kk-1,j+1).*temp1(:));
%                 vminus=sum(v(:,kk-1,j).*ones(p1,1));
%                 penalty=h\beta1*(vplus-vminus);
%                 vxplus=v(2,kk-1,j+1)*2/h+v(3,kk-1,j+1)*2/h*3*(-1)+v(4,kk-1,j+1)*2/h*6;
%                 vxminus=v(2,kk-1,j)*2/h+v(3,kk-1,j)*2/h*3+v(4,kk-1,j)*2/h*6;
%                 vxxplus=v(3,kk-1,j+1)*(2/h)^2*3+v(4,kk-1,j+1)*(2/h)^2*15*(-1);
%                 vxxminus=v(3,kk-1,j)*(2/h)^2*3+v(4,kk-1,j)*(2/h)^2*15;
%                 A2r=Dapp*(penalty+1/2*(vxplus+vxminus)+beta2*h*(vxxplus-vxxminus))-1/2*(vplus+vminus)-Dapp*(1/2*(vplus+vminus)-vminus)*wxminus(i+1);
%                 vplus=sum(v(:,kk-1,j).*temp1(:));
%                 vminus=sum(v(:,kk-1,j-1).*ones(p1,1));
%                 penalty=h\beta1*(vplus-vminus);
%                 vxplus=v(2,kk-1,j)*2/h+v(3,kk-1,j)*2/h*3*(-1)+v(4,kk-1,j)*2/h*6;
%                 vxminus=v(2,kk-1,j-1)*2/h+v(3,kk-1,j-1)*2/h*3+v(4,kk-1,j-1)*2/h*6;
%                 vxxplus=v(3,kk-1,j)*(2/h)^2*3+v(4,kk-1,j)*(2/h)^2*15*(-1);
%                 vxxminus=v(3,kk-1,j-1)*(2/h)^2*3+v(4,kk-1,j-1)*(2/h)^2*15;
%                 A2l=(Dapp*(penalty+1/2*(vxplus+vxminus)+beta2*h*(vxxplus-vxxminus))-1/2*(vplus+vminus))*(-1).^i-Dapp*(1/2*(vplus+vminus)-vplus)*wxplus(i+1);
%             end
%             A1(i+1)=A1r-A1l;
%             A2(i+1)=A2r-A2l;
%         end
%         U(1:2*p1,kk-1,j)=[u(1:p1,kk-1,j);v(1:p1,kk-1,j)];
%         E=[M+E1 F1;F2 M+E2];
%         F=[M+E1+dt*(D_x-Dapp*D_xx) F1;F2 M+E2+dt*(D_x-Dapp*D_xx)];
%         G=[dt*A1';dt*A2'];
%         U(1:2*p1,kk,j)=E\(F*U(1:2*p1,kk-1,j)+G);
%         u(1:p1,kk,j)=U(1:p1,kk,j);
%         v(1:p1,kk,j)=U(p1+1:2*p1,kk,j);
%     end
%     %%% the limiter function down
%     for j=1:J
%         if j==1
%             q=p1;
%             while(q>1)
%                 aa=(2*q-3)*u(q,kk,j);
%                 b=u(q-1,kk,j+1)-u(q-1,kk,j);
%                 %                 c=u(q-1,kk,j)-u(q-1,kk,J);%Periodic boundary conditions
%                 c=u(q-1,kk,j)-umin;%left boundary conditions
%                 d=minmod(aa,b,c);
%                 u(q,kk,j)=(2*q-3)\d;
%                 aa=(2*q-3)*v(q,kk,j);
%                 b=v(q-1,kk,j+1)-v(q-1,kk,j);
%                 %                 c=v(q-1,kk,j)-v(q-1,kk,J);%Periodic boundary conditions
%                 c=v(q-1,kk,j)-vmin;%left boundary conditions
%                 d=minmod(aa,b,c);
%                 v(q,kk,j)=(2*q-3)\d;
%                 q=q-1;
%             end
%         elseif j==J
%             q=p1;
%             while(q>1)
%                 aa=(2*q-3)*u(q,kk,j);
%                 %                 b=u(q-1,kk,1)-u(q-1,kk,j);%Periodic boundary conditions
%                 b=-u(q-1,kk,j);%right boundary conditions
%                 c=u(q-1,kk,j)-u(q-1,kk,j-1);
%                 d=minmod(aa,b,c);
%                 u(q,kk,j)=(2*q-3)\d;
%                 aa=(2*q-3)*v(q,kk,j);
%                 %                 b=v(q-1,kk,1)-v(q-1,kk,j);%Periodic boundary conditions
%                 b=-u(q-1,kk,j);%right boundary conditions
%                 c=v(q-1,kk,j)-v(q-1,kk,j-1);
%                 d=minmod(aa,b,c);
%                 v(q,kk,j)=(2*q-3)\d;
%                 q=q-1;
%             end
%         else
%             q=p1;
%             while(q>1)
%                 aa=(2*q-3)*u(q,kk,j);
%                 b=u(q-1,kk,j+1)-u(q-1,kk,j);
%                 c=u(q-1,kk,j)-u(q-1,kk,j-1);
%                 d=minmod(aa,b,c);
%                 u(q,kk,j)=(2*q-3)\d;
%                 aa=(2*q-3)*v(q,kk,j);
%                 b=v(q-1,kk,j+1)-v(q-1,kk,j);
%                 c=v(q-1,kk,j)-v(q-1,kk,j-1);
%                 d=minmod(aa,b,c);
%                 v(q,kk,j)=(2*q-3)\d;
%                 q=q-1;
%             end
%         end
%     end
%     %%% the limiter function up

    %%% 3-Runge-Kutta method;
    % round one
    for j=1:J
        for i=0:p
            for m=0:p
                sum_loc1=sum(w.*P(i+1,:).*P(m+1,:));
                M(i+1,m+1)=h/2*sum_loc1;
                sum_loc2=sum(w.*D_P(m+1,:).*D_P(i+1,:));
                D_xx(i+1,m+1)=2/h*sum_loc2;
                sum_loc3=sum(w.*P(m+1,:).*D_P(i+1,:));
                D_x(i+1,m+1)=sum_loc3;
                sum_loc1_1=0;
                sum_loc1_2=0;
                sum_loc1_3=0;
                sum_loc1_4=0;
                sum_loc1_5=0;
                sum_loc1_6=0;
                for k=0:p
                    sum_loc1_1=sum_loc1_1+w(k+1).*(2-sin(sum(u(1:p1,kk-1,j).*P(1:p+1,k+1)))).*P(i+1,k+1).*P(m+1,k+1);
                    sum_loc1_2=sum_loc1_2+w(k+1).*0.1.*P(i+1,k+1).*P(m+1,k+1);
                    %             sum_loc1_3=sum_loc1_3+w(k+1).*(-1.5*sum(u(1:p1,kk-1,j).*P(1:p+1,k+1))/(1-sum(u(1:p1,kk-1,j).*P(1:p+1,k+1))+sum(v(1:p1,kk-1,j).*P(1:p+1,k+1)))^2.)*P(i+1,k+1);
                    sum_loc1_4=sum_loc1_4+w(k+1).*0.1.*P(i+1,k+1).*P(m+1,k+1);
                    sum_loc1_5=sum_loc1_5+w(k+1).*(2+cos(sum(v(1:p1,kk-1,j).*P(1:p+1,k+1))))*P(i+1,k+1).*P(m+1,k+1);
                    %             sum_loc1_6=sum_loc1_6+w(k+1).*(-3*sum(v(1:p1,kk-1,j).*P(1:p+1,k+1))/(1-sum(u(1:p1,kk-1,j).*P(1:p+1,k+1))+sum(v(1:p1,kk-1,j).*P(1:p+1,k+1)))^2.)*P(i+1,k+1);
                end
                E1(i+1,m+1)=h/2*sum_loc1_1;
                F1(i+1,m+1)=h/2*sum_loc1_2;
                %                 G1(i+1)=h/2*sum_loc1_3;
                E2(i+1,m+1)=h/2*sum_loc1_4;
                F2(i+1,m+1)=h/2*sum_loc1_5;
                %                 G2(i+1)=h/2*sum_loc1_6;

            end
            if j==1
                uplus=sum(u(:,kk-1,j+1).*temp1(:));
                uminus=sum(u(:,kk-1,j).*ones(p1,1));
                penalty=h\beta1*(uplus-uminus);
                uxplus=u(2,kk-1,j+1)*2/h+u(3,kk-1,j+1)*2/h*3*(-1)+u(4,kk-1,j+1)*2/h*6;
                uxminus=u(2,kk-1,j)*2/h+u(3,kk-1,j)*2/h*3+u(4,kk-1,j)*2/h*6;
                uxxplus=u(3,kk-1,j+1)*(2/h)^2*3+u(4,kk-1,j+1)*(2/h)^2*15*(-1);
                uxxminus=u(3,kk-1,j)*(2/h)^2*3+u(4,kk-1,j)*(2/h)^2*15;
                A1r=Dapp*(penalty+1/2*(uxplus+uxminus)+beta2*h*(uxxplus-uxxminus))-uminus;
                uplus=sum(u(:,kk-1,j).*temp1(:));
                uminus=umin;%left boundary conditions
                penalty=h\beta1*(uplus-uminus);
                uxplus=u(2,kk-1,j)*2/h+u(3,kk-1,j)*2/h*3*(-1)+u(4,kk-1,j)*2/h*6;
                uxminus=0;
                uxxplus=u(3,kk-1,j)*(2/h)^2*3+u(4,kk-1,j)*(2/h)^2*15*(-1);
                uxxminus=0;
                A1l=(Dapp*(penalty+1/2*(uxplus+uxminus)+beta2*h*(uxxplus-uxxminus))-uminus)*(-1).^i;
                vplus=sum(v(:,kk-1,j+1).*temp1(:));
                vminus=sum(v(:,kk-1,j).*ones(p1,1));
                penalty=h\beta1*(vplus-vminus);
                vxplus=v(2,kk-1,j+1)*2/h+v(3,kk-1,j+1)*2/h*3*(-1)+v(4,kk-1,j+1)*2/h*6;
                vxminus=v(2,kk-1,j)*2/h+v(3,kk-1,j)*2/h*3+v(4,kk-1,j)*2/h*6;
                vxxplus=v(3,kk-1,j+1)*(2/h)^2*3+v(4,kk-1,j+1)*(2/h)^2*15*(-1);
                vxxminus=v(3,kk-1,j)*(2/h)^2*3+v(4,kk-1,j)*(2/h)^2*15;
                A2r=Dapp*(penalty+1/2*(vxplus+vxminus)+beta2*h*(vxxplus-vxxminus))-vminus;
                vplus=sum(v(:,kk-1,j).*temp1(:));
                vminus=vmin;%left boundary conditions
                penalty=h\beta1*(vplus-vminus);
                vxplus=v(2,kk-1,j)*2/h+v(3,kk-1,j)*2/h*3*(-1)+v(4,kk-1,j)*2/h*6;
                vxminus=0;
                vxxplus=v(3,kk-1,j)*(2/h)^2*3+v(4,kk-1,j)*(2/h)^2*15*(-1);
                vxxminus=0;
                A2l=(Dapp*(penalty+1/2*(vxplus+vxminus)+beta2*h*(vxxplus-vxxminus))-vminus)*(-1).^i;
            elseif j==J
                uplus=sum(u(:,kk-1,j).*ones(p1,1));%right boundary conditions
                uminus=sum(u(:,kk-1,j).*ones(p1,1));
                penalty=h\beta1*(uplus-uminus);
                uxplus=u(2,kk-1,j)*2/h+u(3,kk-1,j)*2/h*3*(-1)+u(4,kk-1,j)*2/h*6;
                uxminus=u(2,kk-1,j)*2/h+u(3,kk-1,j)*2/h*3+u(4,kk-1,j)*2/h*6;
                uxxplus=u(3,kk-1,j)*(2/h)^2*3+u(4,kk-1,j)*(2/h)^2*15*(-1);
                uxxminus=u(3,kk-1,j)*(2/h)^2*3+u(4,kk-1,j)*(2/h)^2*15;
                A1r=Dapp*(penalty+1/2*(uxplus+uxminus)+beta2*h*(uxxplus-uxxminus))-uminus;
                uplus=sum(u(:,kk-1,j).*temp1(:));
                uminus=sum(u(:,kk-1,j-1).*ones(p1,1));
                penalty=h\beta1*(uplus-uminus);
                uxplus=u(2,kk-1,j)*2/h+u(3,kk-1,j)*2/h*3*(-1)+u(4,kk-1,j)*2/h*6;
                uxminus=u(2,kk-1,j-1)*2/h+u(3,kk-1,j-1)*2/h*3+u(4,kk-1,j-1)*2/h*6;
                uxxplus=u(3,kk-1,j)*(2/h)^2*3+u(4,kk-1,j)*(2/h)^2*15*(-1);
                uxxminus=u(3,kk-1,j-1)*(2/h)^2*3+u(4,kk-1,j-1)*(2/h)^2*15;
                A1l=(Dapp*(penalty+1/2*(uxplus+uxminus)+beta2*h*(uxxplus-uxxminus))-uminus)*(-1).^i;
                vplus=sum(v(:,kk-1,j).*ones(p1,1));%right boundary conditions
                vminus=sum(v(:,kk-1,j).*ones(p1,1));
                penalty=h\beta1*(vplus-vminus);
                vxplus=v(2,kk-1,j)*2/h+v(3,kk-1,j)*2/h*3*(-1)+v(4,kk-1,j)*2/h*6;
                vxminus=v(2,kk-1,j)*2/h+v(3,kk-1,j)*2/h*3+v(4,kk-1,j)*2/h*6;
                vxxplus=v(3,kk-1,j)*(2/h)^2*3+v(4,kk-1,j)*(2/h)^2*15*(-1);
                vxxminus=v(3,kk-1,j)*(2/h)^2*3+v(4,kk-1,j)*(2/h)^2*15;
                A2r=Dapp*(penalty+1/2*(vxplus+vxminus)+beta2*h*(vxxplus-vxxminus))-vminus;
                vplus=sum(v(:,kk-1,j).*temp1(:));
                vminus=sum(v(:,kk-1,j-1).*ones(p1,1));
                penalty=h\beta1*(vplus-vminus);
                vxplus=v(2,kk-1,j)*2/h+v(3,kk-1,j)*2/h*3*(-1)+v(4,kk-1,j)*2/h*6;
                vxminus=v(2,kk-1,j-1)*2/h+v(3,kk-1,j-1)*2/h*3+v(4,kk-1,j-1)*2/h*6;
                vxxplus=v(3,kk-1,j)*(2/h)^2*3+v(4,kk-1,j)*(2/h)^2*15*(-1);
                vxxminus=v(3,kk-1,j-1)*(2/h)^2*3+v(4,kk-1,j-1)*(2/h)^2*15;
                A2l=(Dapp*(penalty+1/2*(vxplus+vxminus)+beta2*h*(vxxplus-vxxminus))-vminus)*(-1).^i;
            else
                uplus=sum(u(:,kk-1,j+1).*temp1(:));
                uminus=sum(u(:,kk-1,j).*ones(p1,1));
                penalty=h\beta1*(uplus-uminus);
                uxplus=u(2,kk-1,j+1)*2/h+u(3,kk-1,j+1)*2/h*3*(-1)+u(4,kk-1,j+1)*2/h*6;
                uxminus=u(2,kk-1,j)*2/h+u(3,kk-1,j)*2/h*3+u(4,kk-1,j)*2/h*6;
                uxxplus=u(3,kk-1,j+1)*(2/h)^2*3+u(4,kk-1,j+1)*(2/h)^2*15*(-1);
                uxxminus=u(3,kk-1,j)*(2/h)^2*3+u(4,kk-1,j)*(2/h)^2*15;
                A1r=Dapp*(penalty+1/2*(uxplus+uxminus)+beta2*h*(uxxplus-uxxminus))-uminus;
                uplus=sum(u(:,kk-1,j).*temp1(:));
                uminus=sum(u(:,kk-1,j-1).*ones(p1,1));
                penalty=h\beta1*(uplus-uminus);
                uxplus=u(2,kk-1,j)*2/h+u(3,kk-1,j)*2/h*3*(-1)+u(4,kk-1,j)*2/h*6;
                uxminus=u(2,kk-1,j-1)*2/h+u(3,kk-1,j-1)*2/h*3+u(4,kk-1,j-1)*2/h*6;
                uxxplus=u(3,kk-1,j)*(2/h)^2*3+u(4,kk-1,j)*(2/h)^2*15*(-1);
                uxxminus=u(3,kk-1,j-1)*(2/h)^2*3+u(4,kk-1,j-1)*(2/h)^2*15;
                A1l=(Dapp*(penalty+1/2*(uxplus+uxminus)+beta2*h*(uxxplus-uxxminus))-uminus)*(-1).^i;
                vplus=sum(v(:,kk-1,j+1).*temp1(:));
                vminus=sum(v(:,kk-1,j).*ones(p1,1));
                penalty=h\beta1*(vplus-vminus);
                vxplus=v(2,kk-1,j+1)*2/h+v(3,kk-1,j+1)*2/h*3*(-1)+v(4,kk-1,j+1)*2/h*6;
                vxminus=v(2,kk-1,j)*2/h+v(3,kk-1,j)*2/h*3+v(4,kk-1,j)*2/h*6;
                vxxplus=v(3,kk-1,j+1)*(2/h)^2*3+v(4,kk-1,j+1)*(2/h)^2*15*(-1);
                vxxminus=v(3,kk-1,j)*(2/h)^2*3+v(4,kk-1,j)*(2/h)^2*15;
                A2r=Dapp*(penalty+1/2*(vxplus+vxminus)+beta2*h*(vxxplus-vxxminus))-vminus;
                vplus=sum(v(:,kk-1,j).*temp1(:));
                vminus=sum(v(:,kk-1,j-1).*ones(p1,1));
                penalty=h\beta1*(vplus-vminus);
                vxplus=v(2,kk-1,j)*2/h+v(3,kk-1,j)*2/h*3*(-1)+v(4,kk-1,j)*2/h*6;
                vxminus=v(2,kk-1,j-1)*2/h+v(3,kk-1,j-1)*2/h*3+v(4,kk-1,j-1)*2/h*6;
                vxxplus=v(3,kk-1,j)*(2/h)^2*3+v(4,kk-1,j)*(2/h)^2*15*(-1);
                vxxminus=v(3,kk-1,j-1)*(2/h)^2*3+v(4,kk-1,j-1)*(2/h)^2*15;
                A2l=(Dapp*(penalty+1/2*(vxplus+vxminus)+beta2*h*(vxxplus-vxxminus))-vminus)*(-1).^i;
            end
            A1(i+1)=A1r-A1l;
            A2(i+1)=A2r-A2l;
        end
        U(1:2*p1,kk-1,j)=[u(1:p1,kk-1,j);v(1:p1,kk-1,j)];
        E=[M+E1 F1;F2 M+E2];
        FF1=zeros(length(F1));
        FF2=zeros(length(F2));
        F=[dt*(D_x-Dapp*D_xx) FF1;FF2 dt*(D_x-Dapp*D_xx)];
        G=[dt*A1';dt*A2'];
        U1(1:2*p1,kk-1,j)=U(1:2*p1,kk-1,j)+E\(F*U(1:2*p1,kk-1,j)+G);
        u1(1:p1,kk-1,j)=U1(1:p1,kk-1,j);
        v1(1:p1,kk-1,j)=U1(p1+1:2*p1,kk-1,j);
    end
    %%% the limiter function down
    for j=1:J
        if j==1
            q=p1;
            while(q>1)
                aa=(2*q-3)*u1(q,kk-1,j);
                b=u1(q-1,kk-1,j+1)-u1(q-1,kk-1,j);
                % c=u1(q-1,kk-1,j)-u1(q-1,kk-1,J);%Periodic boundary conditions
                c=u1(q-1,kk-1,j)-umin;%left boundary conditions
                d=minmod(aa,b,c);
                u1(q,kk-1,j)=(2*q-3)\d;
                aa=(2*q-3)*v1(q,kk-1,j);
                b=v1(q-1,kk-1,j+1)-v1(q-1,kk-1,j);
                % c=v1(q-1,kk-1,j)-v1(q-1,kk-1,J);%Periodic boundary conditions
                c=v1(q-1,kk-1,j)-vmin;%left boundary conditions
                d=minmod(aa,b,c);
                v1(q,kk-1,j)=(2*q-3)\d;
                q=q-1;
            end
        elseif j==J
            q=p1;
            while(q>1)
                aa=(2*q-3)*u1(q,kk-1,j);
                % b=u1(q-1,kk-1,1)-u1(q-1,kk-1,j);%Periodic boundary conditions
                b=-u1(q-1,kk-1,j);%right boundary conditions
                c=u1(q-1,kk-1,j)-u1(q-1,kk-1,j-1);
                d=minmod(aa,b,c);
                u1(q,kk-1,j)=(2*q-3)\d;
                aa=(2*q-3)*v1(q,kk-1,j);
                % b=v1(q-1,kk-1,1)-v1(q-1,kk-1,j);%Periodic boundary conditions
                b=-u1(q-1,kk-1,j);%right boundary conditions
                c=v1(q-1,kk-1,j)-v1(q-1,kk-1,j-1);
                d=minmod(aa,b,c);
                v1(q,kk-1,j)=(2*q-3)\d;
                q=q-1;
            end
        else
            q=p1;
            while(q>1)
                aa=(2*q-3)*u1(q,kk-1,j);
                b=u1(q-1,kk-1,j+1)-u1(q-1,kk-1,j);
                c=u1(q-1,kk-1,j)-u1(q-1,kk-1,j-1);
                d=minmod(aa,b,c);
                u1(q,kk-1,j)=(2*q-3)\d;
                aa=(2*q-3)*v1(q,kk-1,j);
                b=v1(q-1,kk-1,j+1)-v1(q-1,kk-1,j);
                c=v1(q-1,kk-1,j)-v1(q-1,kk-1,j-1);
                d=minmod(aa,b,c);
                v1(q,kk-1,j)=(2*q-3)\d;
                q=q-1;
            end
        end
    end
    %%% the limiter function up

    % round two
    for j=1:J
        for i=0:p
            for m=0:p
                sum_loc1=sum(w.*P(i+1,:).*P(m+1,:));
                M(i+1,m+1)=h/2*sum_loc1;
                sum_loc2=sum(w.*D_P(m+1,:).*D_P(i+1,:));
                D_xx(i+1,m+1)=2/h*sum_loc2;
                sum_loc3=sum(w.*P(m+1,:).*D_P(i+1,:));
                D_x(i+1,m+1)=sum_loc3;
                sum_loc1_1=0;
                sum_loc1_2=0;
                sum_loc1_3=0;
                sum_loc1_4=0;
                sum_loc1_5=0;
                sum_loc1_6=0;
                for k=0:p
                    sum_loc1_1=sum_loc1_1+w(k+1).*(2-sin(sum(u1(1:p1,kk-1,j).*P(1:p+1,k+1)))).*P(i+1,k+1).*P(m+1,k+1);
                    sum_loc1_2=sum_loc1_2+w(k+1).*0.1.*P(i+1,k+1).*P(m+1,k+1);
                    %             sum_loc1_3=sum_loc1_3+w(k+1).*(-1.5*sum(u(1:p1,kk-1,j).*P(1:p+1,k+1))/(1-sum(u(1:p1,kk-1,j).*P(1:p+1,k+1))+sum(v(1:p1,kk-1,j).*P(1:p+1,k+1)))^2.)*P(i+1,k+1);
                    sum_loc1_4=sum_loc1_4+w(k+1).*0.1.*P(i+1,k+1).*P(m+1,k+1);
                    sum_loc1_5=sum_loc1_5+w(k+1).*(2+cos(sum(v1(1:p1,kk-1,j).*P(1:p+1,k+1))))*P(i+1,k+1).*P(m+1,k+1);
                    %             sum_loc1_6=sum_loc1_6+w(k+1).*(-3*sum(v(1:p1,kk-1,j).*P(1:p+1,k+1))/(1-sum(u(1:p1,kk-1,j).*P(1:p+1,k+1))+sum(v(1:p1,kk-1,j).*P(1:p+1,k+1)))^2.)*P(i+1,k+1);
                end
                E1(i+1,m+1)=h/2*sum_loc1_1;
                F1(i+1,m+1)=h/2*sum_loc1_2;
                %                 G1(i+1)=h/2*sum_loc1_3;
                E2(i+1,m+1)=h/2*sum_loc1_4;
                F2(i+1,m+1)=h/2*sum_loc1_5;
                %                 G2(i+1)=h/2*sum_loc1_6;
            end

            if j==1
                uplus=sum(u1(:,kk-1,j+1).*temp1(:));
                uminus=sum(u1(:,kk-1,j).*ones(p1,1));
                penalty=h\beta1*(uplus-uminus);
                uxplus=u1(2,kk-1,j+1)*2/h+u1(3,kk-1,j+1)*2/h*3*(-1)+u1(4,kk-1,j+1)*2/h*6;
                uxminus=u1(2,kk-1,j)*2/h+u1(3,kk-1,j)*2/h*3+u1(4,kk-1,j)*2/h*6;
                uxxplus=u1(3,kk-1,j+1)*(2/h)^2*3+u1(4,kk-1,j+1)*(2/h)^2*15*(-1);
                uxxminus=u1(3,kk-1,j)*(2/h)^2*3+u1(4,kk-1,j)*(2/h)^2*15;
                A1r=Dapp*(penalty+1/2*(uxplus+uxminus)+beta2*h*(uxxplus-uxxminus))-uminus;
                uplus=sum(u1(:,kk-1,j).*temp1(:));
                uminus=umin;%left boundary conditions
                penalty=h\beta1*(uplus-uminus);
                uxplus=u1(2,kk-1,j)*2/h+u1(3,kk-1,j)*2/h*3*(-1)+u1(4,kk-1,j)*2/h*6;
                uxminus=0;
                uxxplus=u1(3,kk-1,j)*(2/h)^2*3+u1(4,kk-1,j)*(2/h)^2*15*(-1);
                uxxminus=0;
                A1l=(Dapp*(penalty+1/2*(uxplus+uxminus)+beta2*h*(uxxplus-uxxminus))-uminus)*(-1).^i;
                vplus=sum(v1(:,kk-1,j+1).*temp1(:));
                vminus=sum(v1(:,kk-1,j).*ones(p1,1));
                penalty=h\beta1*(vplus-vminus);
                vxplus=v1(2,kk-1,j+1)*2/h+v1(3,kk-1,j+1)*2/h*3*(-1)+v1(4,kk-1,j+1)*2/h*6;
                vxminus=v1(2,kk-1,j)*2/h+v1(3,kk-1,j)*2/h*3+v1(4,kk-1,j)*2/h*6;
                vxxplus=v1(3,kk-1,j+1)*(2/h)^2*3+v1(4,kk-1,j+1)*(2/h)^2*15*(-1);
                vxxminus=v1(3,kk-1,j)*(2/h)^2*3+v1(4,kk-1,j)*(2/h)^2*15;
                A2r=Dapp*(penalty+1/2*(vxplus+vxminus)+beta2*h*(vxxplus-vxxminus))-vminus;
                vplus=sum(v1(:,kk-1,j).*temp1(:));
                vminus=vmin;%left boundary conditions
                penalty=h\beta1*(vplus-vminus);
                vxplus=v1(2,kk-1,j)*2/h+v1(3,kk-1,j)*2/h*3*(-1)+v1(4,kk-1,j)*2/h*6;
                vxminus=0;
                vxxplus=v1(3,kk-1,j)*(2/h)^2*3+v1(4,kk-1,j)*(2/h)^2*15*(-1);
                vxxminus=0;
                A2l=(Dapp*(penalty+1/2*(vxplus+vxminus)+beta2*h*(vxxplus-vxxminus))-vminus)*(-1).^i;
            elseif j==J
                uplus=sum(u1(:,kk-1,j).*ones(p1,1));%right boundary conditions
                uminus=sum(u1(:,kk-1,j).*ones(p1,1));
                penalty=h\beta1*(uplus-uminus);
                uxplus=u1(2,kk-1,j)*2/h+u1(3,kk-1,j)*2/h*3*(-1)+u1(4,kk-1,j)*2/h*6;
                uxminus=u1(2,kk-1,j)*2/h+u1(3,kk-1,j)*2/h*3+u1(4,kk-1,j)*2/h*6;
                uxxplus=u1(3,kk-1,j)*(2/h)^2*3+u1(4,kk-1,j)*(2/h)^2*15*(-1);
                uxxminus=u1(3,kk-1,j)*(2/h)^2*3+u1(4,kk-1,j)*(2/h)^2*15;
                A1r=Dapp*(penalty+1/2*(uxplus+uxminus)+beta2*h*(uxxplus-uxxminus))-uminus;
                uplus=sum(u1(:,kk-1,j).*temp1(:));
                uminus=sum(u1(:,kk-1,j-1).*ones(p1,1));
                penalty=h\beta1*(uplus-uminus);
                uxplus=u1(2,kk-1,j)*2/h+u1(3,kk-1,j)*2/h*3*(-1)+u1(4,kk-1,j)*2/h*6;
                uxminus=u1(2,kk-1,j-1)*2/h+u1(3,kk-1,j-1)*2/h*3+u1(4,kk-1,j-1)*2/h*6;
                uxxplus=u1(3,kk-1,j)*(2/h)^2*3+u1(4,kk-1,j)*(2/h)^2*15*(-1);
                uxxminus=u1(3,kk-1,j-1)*(2/h)^2*3+u1(4,kk-1,j-1)*(2/h)^2*15;
                A1l=(Dapp*(penalty+1/2*(uxplus+uxminus)+beta2*h*(uxxplus-uxxminus))-uminus)*(-1).^i;
                vplus=sum(v1(:,kk-1,j).*ones(p1,1));%right boundary conditions
                vminus=sum(v1(:,kk-1,j).*ones(p1,1));
                penalty=h\beta1*(vplus-vminus);
                vxplus=v1(2,kk-1,j)*2/h+v1(3,kk-1,j)*2/h*3*(-1)+v1(4,kk-1,j)*2/h*6;
                vxminus=v1(2,kk-1,j)*2/h+v1(3,kk-1,j)*2/h*3+v1(4,kk-1,j)*2/h*6;
                vxxplus=v1(3,kk-1,j)*(2/h)^2*3+v1(4,kk-1,j)*(2/h)^2*15*(-1);
                vxxminus=v1(3,kk-1,j)*(2/h)^2*3+v1(4,kk-1,j)*(2/h)^2*15;
                A2r=Dapp*(penalty+1/2*(vxplus+vxminus)+beta2*h*(vxxplus-vxxminus))-vminus;
                vplus=sum(v1(:,kk-1,j).*temp1(:));
                vminus=sum(v1(:,kk-1,j-1).*ones(p1,1));
                penalty=h\beta1*(vplus-vminus);
                vxplus=v1(2,kk-1,j)*2/h+v1(3,kk-1,j)*2/h*3*(-1)+v1(4,kk-1,j)*2/h*6;
                vxminus=v1(2,kk-1,j-1)*2/h+v1(3,kk-1,j-1)*2/h*3+v1(4,kk-1,j-1)*2/h*6;
                vxxplus=v1(3,kk-1,j)*(2/h)^2*3+v1(4,kk-1,j)*(2/h)^2*15*(-1);
                vxxminus=v1(3,kk-1,j-1)*(2/h)^2*3+v1(4,kk-1,j-1)*(2/h)^2*15;
                A2l=(Dapp*(penalty+1/2*(vxplus+vxminus)+beta2*h*(vxxplus-vxxminus))-vminus)*(-1).^i;
            else
                uplus=sum(u1(:,kk-1,j+1).*temp1(:));
                uminus=sum(u1(:,kk-1,j).*ones(p1,1));
                penalty=h\beta1*(uplus-uminus);
                uxplus=u1(2,kk-1,j+1)*2/h+u1(3,kk-1,j+1)*2/h*3*(-1)+u1(4,kk-1,j+1)*2/h*6;
                uxminus=u1(2,kk-1,j)*2/h+u1(3,kk-1,j)*2/h*3+u1(4,kk-1,j)*2/h*6;
                uxxplus=u1(3,kk-1,j+1)*(2/h)^2*3+u1(4,kk-1,j+1)*(2/h)^2*15*(-1);
                uxxminus=u1(3,kk-1,j)*(2/h)^2*3+u1(4,kk-1,j)*(2/h)^2*15;
                A1r=Dapp*(penalty+1/2*(uxplus+uxminus)+beta2*h*(uxxplus-uxxminus))-uminus;
                uplus=sum(u1(:,kk-1,j).*temp1(:));
                uminus=sum(u1(:,kk-1,j-1).*ones(p1,1));
                penalty=h\beta1*(uplus-uminus);
                uxplus=u1(2,kk-1,j)*2/h+u1(3,kk-1,j)*2/h*3*(-1)+u1(4,kk-1,j)*2/h*6;
                uxminus=u1(2,kk-1,j-1)*2/h+u1(3,kk-1,j-1)*2/h*3+u1(4,kk-1,j-1)*2/h*6;
                uxxplus=u1(3,kk-1,j)*(2/h)^2*3+u1(4,kk-1,j)*(2/h)^2*15*(-1);
                uxxminus=u1(3,kk-1,j-1)*(2/h)^2*3+u1(4,kk-1,j-1)*(2/h)^2*15;
                A1l=(Dapp*(penalty+1/2*(uxplus+uxminus)+beta2*h*(uxxplus-uxxminus))-uminus)*(-1).^i;
                vplus=sum(v1(:,kk-1,j+1).*temp1(:));
                vminus=sum(v1(:,kk-1,j).*ones(p1,1));
                penalty=h\beta1*(vplus-vminus);
                vxplus=v1(2,kk-1,j+1)*2/h+v1(3,kk-1,j+1)*2/h*3*(-1)+v1(4,kk-1,j+1)*2/h*6;
                vxminus=v1(2,kk-1,j)*2/h+v1(3,kk-1,j)*2/h*3+v1(4,kk-1,j)*2/h*6;
                vxxplus=v1(3,kk-1,j+1)*(2/h)^2*3+v1(4,kk-1,j+1)*(2/h)^2*15*(-1);
                vxxminus=v1(3,kk-1,j)*(2/h)^2*3+v1(4,kk-1,j)*(2/h)^2*15;
                A2r=Dapp*(penalty+1/2*(vxplus+vxminus)+beta2*h*(vxxplus-vxxminus))-vminus;
                vplus=sum(v1(:,kk-1,j).*temp1(:));
                vminus=sum(v1(:,kk-1,j-1).*ones(p1,1));
                penalty=h\beta1*(vplus-vminus);
                vxplus=v1(2,kk-1,j)*2/h+v1(3,kk-1,j)*2/h*3*(-1)+v1(4,kk-1,j)*2/h*6;
                vxminus=v1(2,kk-1,j-1)*2/h+v1(3,kk-1,j-1)*2/h*3+v1(4,kk-1,j-1)*2/h*6;
                vxxplus=v1(3,kk-1,j)*(2/h)^2*3+v1(4,kk-1,j)*(2/h)^2*15*(-1);
                vxxminus=v1(3,kk-1,j-1)*(2/h)^2*3+v1(4,kk-1,j-1)*(2/h)^2*15;
                A2l=(Dapp*(penalty+1/2*(vxplus+vxminus)+beta2*h*(vxxplus-vxxminus))-vminus)*(-1).^i;
            end
            A1(i+1)=A1r-A1l;
            A2(i+1)=A2r-A2l;
        end
        U1(1:2*p1,kk-1,j)=[u1(1:p1,kk-1,j);v1(1:p1,kk-1,j)];
        E=[M+E1 F1;F2 M+E2];
        FF1=zeros(length(F1));
        FF2=zeros(length(F2));
        F=[dt*(D_x-Dapp*D_xx) FF1;FF2 dt*(D_x-Dapp*D_xx)];
        G=[dt*A1';dt*A2'];
        U2(1:2*p1,kk-1,j)=3/4*U(1:2*p1,kk-1,j)+1/4*(U1(1:2*p1,kk-1,j)+E\(F*U1(1:2*p1,kk-1,j)+G));
        u2(1:p1,kk-1,j)=U2(1:p1,kk-1,j);
        v2(1:p1,kk-1,j)=U2(p1+1:2*p1,kk-1,j);
    end
    %%% the limiter function down
    for j=1:J
        if j==1
            q=p1;
            while(q>1)
                aa=(2*q-3)*u2(q,kk-1,j);
                b=u2(q-1,kk-1,j+1)-u2(q-1,kk-1,j);
                % c=u2(q-1,kk-1,j)-u2(q-1,kk-1,J);%Periodic boundary conditions
                c=u2(q-1,kk-1,j)-umin;%left boundary conditions
                d=minmod(aa,b,c);
                u2(q,kk-1,j)=(2*q-3)\d;
                aa=(2*q-3)*v2(q,kk-1,j);
                b=v2(q-1,kk-1,j+1)-v2(q-1,kk-1,j);
                % c=v2(q-1,kk-1,j)-v2(q-1,kk-1,J);%Periodic boundary conditions
                c=v2(q-1,kk-1,j)-vmin;%left boundary conditions
                d=minmod(aa,b,c);
                v2(q,kk-1,j)=(2*q-3)\d;
                q=q-1;
            end
        elseif j==J
            q=p1;
            while(q>1)
                aa=(2*q-3)*u2(q,kk-1,j);
                % b=u2(q-1,kk-1,1)-u2(q-1,kk-1,j);%Periodic boundary conditions
                b=-u2(q-1,kk-1,j);%right boundary conditions
                c=u2(q-1,kk-1,j)-u2(q-1,kk-1,j-1);
                d=minmod(aa,b,c);
                u2(q,kk-1,j)=(2*q-3)\d;
                aa=(2*q-3)*v2(q,kk-1,j);
                % b=v2(q-1,kk-1,1)-v2(q-1,kk-1,j);%Periodic boundary conditions
                b=-u2(q-1,kk-1,j);%right boundary conditions
                c=v2(q-1,kk-1,j)-v2(q-1,kk-1,j-1);
                d=minmod(aa,b,c);
                v2(q,kk-1,j)=(2*q-3)\d;
                q=q-1;
            end
        else
            q=p1;
            while(q>1)
                aa=(2*q-3)*u2(q,kk-1,j);
                b=u2(q-1,kk-1,j+1)-u2(q-1,kk-1,j);
                c=u2(q-1,kk-1,j)-u2(q-1,kk-1,j-1);
                d=minmod(aa,b,c);
                u2(q,kk-1,j)=(2*q-3)\d;
                aa=(2*q-3)*v2(q,kk-1,j);
                b=v2(q-1,kk-1,j+1)-v2(q-1,kk-1,j);
                c=v2(q-1,kk-1,j)-v2(q-1,kk-1,j-1);
                d=minmod(aa,b,c);
                v2(q,kk-1,j)=(2*q-3)\d;
                q=q-1;
            end
        end
    end
    %%% the limiter function up
    % round three
    for j=1:J
        for i=0:p
            for m=0:p
                sum_loc1=sum(w.*P(i+1,:).*P(m+1,:));
                M(i+1,m+1)=h/2*sum_loc1;
                sum_loc2=sum(w.*D_P(m+1,:).*D_P(i+1,:));
                D_xx(i+1,m+1)=2/h*sum_loc2;
                sum_loc3=sum(w.*P(m+1,:).*D_P(i+1,:));
                D_x(i+1,m+1)=sum_loc3;
                sum_loc1_1=0;
                sum_loc1_2=0;
                sum_loc1_3=0;
                sum_loc1_4=0;
                sum_loc1_5=0;
                sum_loc1_6=0;
                for k=0:p
                    sum_loc1_1=sum_loc1_1+w(k+1).*(2-sin(sum(u2(1:p1,kk-1,j).*P(1:p+1,k+1)))).*P(i+1,k+1).*P(m+1,k+1);
                    sum_loc1_2=sum_loc1_2+w(k+1).*0.1.*P(i+1,k+1).*P(m+1,k+1);
                    %             sum_loc1_3=sum_loc1_3+w(k+1).*(-1.5*sum(u(1:p1,kk-1,j).*P(1:p+1,k+1))/(1-sum(u(1:p1,kk-1,j).*P(1:p+1,k+1))+sum(v(1:p1,kk-1,j).*P(1:p+1,k+1)))^2.)*P(i+1,k+1);
                    sum_loc1_4=sum_loc1_4+w(k+1).*0.1.*P(i+1,k+1).*P(m+1,k+1);
                    sum_loc1_5=sum_loc1_5+w(k+1).*(2+cos(sum(v2(1:p1,kk-1,j).*P(1:p+1,k+1))))*P(i+1,k+1).*P(m+1,k+1);
                    %             sum_loc1_6=sum_loc1_6+w(k+1).*(-3*sum(v(1:p1,kk-1,j).*P(1:p+1,k+1))/(1-sum(u(1:p1,kk-1,j).*P(1:p+1,k+1))+sum(v(1:p1,kk-1,j).*P(1:p+1,k+1)))^2.)*P(i+1,k+1);
                end
                E1(i+1,m+1)=h/2*sum_loc1_1;
                F1(i+1,m+1)=h/2*sum_loc1_2;
                %                 G1(i+1)=h/2*sum_loc1_3;
                E2(i+1,m+1)=h/2*sum_loc1_4;
                F2(i+1,m+1)=h/2*sum_loc1_5;
                %                 G2(i+1)=h/2*sum_loc1_6;
            end

            if j==1
                uplus=sum(u2(:,kk-1,j+1).*temp1(:));
                uminus=sum(u2(:,kk-1,j).*ones(p1,1));
                penalty=h\beta1*(uplus-uminus);
                uxplus=u2(2,kk-1,j+1)*2/h+u2(3,kk-1,j+1)*2/h*3*(-1)+u2(4,kk-1,j+1)*2/h*6;
                uxminus=u2(2,kk-1,j)*2/h+u2(3,kk-1,j)*2/h*3+u2(4,kk-1,j)*2/h*6;
                uxxplus=u2(3,kk-1,j+1)*(2/h)^2*3+u2(4,kk-1,j+1)*(2/h)^2*15*(-1);
                uxxminus=u2(3,kk-1,j)*(2/h)^2*3+u2(4,kk-1,j)*(2/h)^2*15;
                A1r=Dapp*(penalty+1/2*(uxplus+uxminus)+beta2*h*(uxxplus-uxxminus))-uminus;
                uplus=sum(u2(:,kk-1,j).*temp1(:));
                uminus=umin;%left boundary conditions
                penalty=h\beta1*(uplus-uminus);
                uxplus=u2(2,kk-1,j)*2/h+u2(3,kk-1,j)*2/h*3*(-1)+u2(4,kk-1,j)*2/h*6;
                uxminus=0;
                uxxplus=u2(3,kk-1,j)*(2/h)^2*3+u2(4,kk-1,j)*(2/h)^2*15*(-1);
                uxxminus=0;
                A1l=(Dapp*(penalty+1/2*(uxplus+uxminus)+beta2*h*(uxxplus-uxxminus))-uminus)*(-1).^i;
                vplus=sum(v2(:,kk-1,j+1).*temp1(:));
                vminus=sum(v2(:,kk-1,j).*ones(p1,1));
                penalty=h\beta1*(vplus-vminus);
                vxplus=v2(2,kk-1,j+1)*2/h+v2(3,kk-1,j+1)*2/h*3*(-1)+v2(4,kk-1,j+1)*2/h*6;
                vxminus=v2(2,kk-1,j)*2/h+v2(3,kk-1,j)*2/h*3+v2(4,kk-1,j)*2/h*6;
                vxxplus=v2(3,kk-1,j+1)*(2/h)^2*3+v2(4,kk-1,j+1)*(2/h)^2*15*(-1);
                vxxminus=v2(3,kk-1,j)*(2/h)^2*3+v2(4,kk-1,j)*(2/h)^2*15;
                A2r=Dapp*(penalty+1/2*(vxplus+vxminus)+beta2*h*(vxxplus-vxxminus))-vminus;
                vplus=sum(v2(:,kk-1,j).*temp1(:));
                vminus=vmin;%left boundary conditions
                penalty=h\beta1*(vplus-vminus);
                vxplus=v2(2,kk-1,j)*2/h+v2(3,kk-1,j)*2/h*3*(-1)+v2(4,kk-1,j)*2/h*6;
                vxminus=0;
                vxxplus=v2(3,kk-1,j)*(2/h)^2*3+v2(4,kk-1,j)*(2/h)^2*15*(-1);
                vxxminus=0;
                A2l=(Dapp*(penalty+1/2*(vxplus+vxminus)+beta2*h*(vxxplus-vxxminus))-vminus)*(-1).^i;
            elseif j==J
                uplus=sum(u2(:,kk-1,j).*ones(p1,1));%right boundary conditions
                uminus=sum(u2(:,kk-1,j).*ones(p1,1));
                penalty=h\beta1*(uplus-uminus);
                uxplus=u2(2,kk-1,j)*2/h+u2(3,kk-1,j)*2/h*3*(-1)+u2(4,kk-1,j)*2/h*6;
                uxminus=u2(2,kk-1,j)*2/h+u2(3,kk-1,j)*2/h*3+u2(4,kk-1,j)*2/h*6;
                uxxplus=u2(3,kk-1,j)*(2/h)^2*3+u2(4,kk-1,j)*(2/h)^2*15*(-1);
                uxxminus=u2(3,kk-1,j)*(2/h)^2*3+u2(4,kk-1,j)*(2/h)^2*15;
                A1r=Dapp*(penalty+1/2*(uxplus+uxminus)+beta2*h*(uxxplus-uxxminus))-uminus;
                uplus=sum(u2(:,kk-1,j).*temp1(:));
                uminus=sum(u2(:,kk-1,j-1).*ones(p1,1));
                penalty=h\beta1*(uplus-uminus);
                uxplus=u2(2,kk-1,j)*2/h+u2(3,kk-1,j)*2/h*3*(-1)+u2(4,kk-1,j)*2/h*6;
                uxminus=u2(2,kk-1,j-1)*2/h+u2(3,kk-1,j-1)*2/h*3+u2(4,kk-1,j-1)*2/h*6;
                uxxplus=u2(3,kk-1,j)*(2/h)^2*3+u2(4,kk-1,j)*(2/h)^2*15*(-1);
                uxxminus=u2(3,kk-1,j-1)*(2/h)^2*3+u2(4,kk-1,j-1)*(2/h)^2*15;
                A1l=(Dapp*(penalty+1/2*(uxplus+uxminus)+beta2*h*(uxxplus-uxxminus))-uminus)*(-1).^i;
                vplus=sum(v2(:,kk-1,j).*ones(p1,1));%right boundary conditions
                vminus=sum(v2(:,kk-1,j).*ones(p1,1));
                penalty=h\beta1*(vplus-vminus);
                vxplus=v2(2,kk-1,j)*2/h+v2(3,kk-1,j)*2/h*3*(-1)+v2(4,kk-1,j)*2/h*6;
                vxminus=v2(2,kk-1,j)*2/h+v2(3,kk-1,j)*2/h*3+v2(4,kk-1,j)*2/h*6;
                vxxplus=v2(3,kk-1,j)*(2/h)^2*3+v2(4,kk-1,j)*(2/h)^2*15*(-1);
                vxxminus=v2(3,kk-1,j)*(2/h)^2*3+v2(4,kk-1,j)*(2/h)^2*15;
                A2r=Dapp*(penalty+1/2*(vxplus+vxminus)+beta2*h*(vxxplus-vxxminus))-vminus;
                vplus=sum(v2(:,kk-1,j).*temp1(:));
                vminus=sum(v2(:,kk-1,j-1).*ones(p1,1));
                penalty=h\beta1*(vplus-vminus);
                vxplus=v2(2,kk-1,j)*2/h+v2(3,kk-1,j)*2/h*3*(-1)+v2(4,kk-1,j)*2/h*6;
                vxminus=v2(2,kk-1,j-1)*2/h+v2(3,kk-1,j-1)*2/h*3+v2(4,kk-1,j-1)*2/h*6;
                vxxplus=v2(3,kk-1,j)*(2/h)^2*3+v2(4,kk-1,j)*(2/h)^2*15*(-1);
                vxxminus=v2(3,kk-1,j-1)*(2/h)^2*3+v2(4,kk-1,j-1)*(2/h)^2*15;
                A2l=(Dapp*(penalty+1/2*(vxplus+vxminus)+beta2*h*(vxxplus-vxxminus))-vminus)*(-1).^i;
            else
                uplus=sum(u2(:,kk-1,j+1).*temp1(:));
                uminus=sum(u2(:,kk-1,j).*ones(p1,1));
                penalty=h\beta1*(uplus-uminus);
                uxplus=u2(2,kk-1,j+1)*2/h+u2(3,kk-1,j+1)*2/h*3*(-1)+u2(4,kk-1,j+1)*2/h*6;
                uxminus=u2(2,kk-1,j)*2/h+u2(3,kk-1,j)*2/h*3+u2(4,kk-1,j)*2/h*6;
                uxxplus=u2(3,kk-1,j+1)*(2/h)^2*3+u2(4,kk-1,j+1)*(2/h)^2*15*(-1);
                uxxminus=u2(3,kk-1,j)*(2/h)^2*3+u2(4,kk-1,j)*(2/h)^2*15;
                A1r=Dapp*(penalty+1/2*(uxplus+uxminus)+beta2*h*(uxxplus-uxxminus))-uminus;
                uplus=sum(u2(:,kk-1,j).*temp1(:));
                uminus=sum(u2(:,kk-1,j-1).*ones(p1,1));
                penalty=h\beta1*(uplus-uminus);
                uxplus=u2(2,kk-1,j)*2/h+u2(3,kk-1,j)*2/h*3*(-1)+u2(4,kk-1,j)*2/h*6;
                uxminus=u2(2,kk-1,j-1)*2/h+u2(3,kk-1,j-1)*2/h*3+u2(4,kk-1,j-1)*2/h*6;
                uxxplus=u2(3,kk-1,j)*(2/h)^2*3+u2(4,kk-1,j)*(2/h)^2*15*(-1);
                uxxminus=u2(3,kk-1,j-1)*(2/h)^2*3+u2(4,kk-1,j-1)*(2/h)^2*15;
                A1l=(Dapp*(penalty+1/2*(uxplus+uxminus)+beta2*h*(uxxplus-uxxminus))-uminus)*(-1).^i;
                vplus=sum(v2(:,kk-1,j+1).*temp1(:));
                vminus=sum(v2(:,kk-1,j).*ones(p1,1));
                penalty=h\beta1*(vplus-vminus);
                vxplus=v2(2,kk-1,j+1)*2/h+v2(3,kk-1,j+1)*2/h*3*(-1)+v2(4,kk-1,j+1)*2/h*6;
                vxminus=v2(2,kk-1,j)*2/h+v2(3,kk-1,j)*2/h*3+v2(4,kk-1,j)*2/h*6;
                vxxplus=v2(3,kk-1,j+1)*(2/h)^2*3+v2(4,kk-1,j+1)*(2/h)^2*15*(-1);
                vxxminus=v2(3,kk-1,j)*(2/h)^2*3+v2(4,kk-1,j)*(2/h)^2*15;
                A2r=Dapp*(penalty+1/2*(vxplus+vxminus)+beta2*h*(vxxplus-vxxminus))-vminus;
                vplus=sum(v2(:,kk-1,j).*temp1(:));
                vminus=sum(v2(:,kk-1,j-1).*ones(p1,1));
                penalty=h\beta1*(vplus-vminus);
                vxplus=v2(2,kk-1,j)*2/h+v2(3,kk-1,j)*2/h*3*(-1)+v2(4,kk-1,j)*2/h*6;
                vxminus=v2(2,kk-1,j-1)*2/h+v2(3,kk-1,j-1)*2/h*3+v2(4,kk-1,j-1)*2/h*6;
                vxxplus=v2(3,kk-1,j)*(2/h)^2*3+v2(4,kk-1,j)*(2/h)^2*15*(-1);
                vxxminus=v2(3,kk-1,j-1)*(2/h)^2*3+v2(4,kk-1,j-1)*(2/h)^2*15;
                A2l=(Dapp*(penalty+1/2*(vxplus+vxminus)+beta2*h*(vxxplus-vxxminus))-vminus)*(-1).^i;
            end
            A1(i+1)=A1r-A1l;
            A2(i+1)=A2r-A2l;
        end
        U2(1:2*p1,kk-1,j)=[u2(1:p1,kk-1,j);v2(1:p1,kk-1,j)];
        E=[M+E1 F1;F2 M+E2];
        FF1=zeros(length(F1));
        FF2=zeros(length(F2));
        F=[dt*(D_x-Dapp*D_xx) FF1;FF2 dt*(D_x-Dapp*D_xx)];
        G=[dt*A1';dt*A2'];
        U(1:2*p1,kk,j)=1/3*U(1:2*p1,kk-1,j)+2/3*(U2(1:2*p1,kk-1,j)+E\(F*U2(1:2*p1,kk-1,j)+G));
        u(1:p1,kk,j)=U(1:p1,kk,j);
        v(1:p1,kk,j)=U(p1+1:2*p1,kk,j);
    end
    %%% the limiter function down
    for j=1:J
        if j==1
            q=p1;
            while(q>1)
                aa=(2*q-3)*u(q,kk,j);
                b=u(q-1,kk,j+1)-u(q-1,kk,j);
                % c=u(q-1,kk,j)-u(q-1,kk,J);%Periodic boundary conditions
                c=u(q-1,kk,j)-umin;%left boundary conditions
                d=minmod(aa,b,c);
                u(q,kk,j)=(2*q-3)\d;
                aa=(2*q-3)*v(q,kk,j);
                b=v(q-1,kk,j+1)-v(q-1,kk,j);
                % c=v(q-1,kk,j)-v(q-1,kk,J);%Periodic boundary conditions
                c=v(q-1,kk,j)-vmin;%left boundary conditions
                d=minmod(aa,b,c);
                v(q,kk,j)=(2*q-3)\d;
                q=q-1;
            end
        elseif j==J
            q=p1;
            while(q>1)
                aa=(2*q-3)*u(q,kk,j);
                % b=u(q-1,kk,1)-u(q-1,kk,j);%Periodic boundary conditions
                b=-u(q-1,kk,j);%right boundary conditions
                c=u(q-1,kk,j)-u(q-1,kk,j-1);
                d=minmod(aa,b,c);
                u(q,kk,j)=(2*q-3)\d;
                aa=(2*q-3)*v(q,kk,j);
                % b=v(q-1,kk,1)-v(q-1,kk,j);%Periodic boundary conditions
                b=-u(q-1,kk,j);%right boundary conditions
                c=v(q-1,kk,j)-v(q-1,kk,j-1);
                d=minmod(aa,b,c);
                v(q,kk,j)=(2*q-3)\d;
                q=q-1;
            end
        else
            q=p1;
            while(q>1)
                aa=(2*q-3)*u(q,kk,j);
                b=u(q-1,kk,j+1)-u(q-1,kk,j);
                c=u(q-1,kk,j)-u(q-1,kk,j-1);
                d=minmod(aa,b,c);
                u(q,kk,j)=(2*q-3)\d;
                aa=(2*q-3)*v(q,kk,j);
                b=v(q-1,kk,j+1)-v(q-1,kk,j);
                c=v(q-1,kk,j)-v(q-1,kk,j-1);
                d=minmod(aa,b,c);
                v(q,kk,j)=(2*q-3)\d;
                q=q-1;
            end
        end
    end
    %%% the limiter function up

    time=time+dt;
end
%% return the value of function u
for k=1:kk
    for j=1:J
        unum(k,1+(j-1)*p1:j*p1)=PP*u(:,k,j);
        vnum(k,1+(j-1)*p1:j*p1)=PP*v(:,k,j);
    end
end
%% x-coordinate
for j=1:J
    xx(1+(j-1)*p1:j*p1)=x_loc(j,2:p1+1);
end
%% plot the result
%dt=final_t/Nt;
dt=deltat;
time=-dt;
kk=0;
while(time<final_t)
    if time+dt>final_t
        dt=final_t-time;
    end
    kk=kk+1;
    time=time+dt;
    unum1(kk)=unum(kk,end);
    vnum1(kk)=vnum(kk,end);
    plot(xx,unum(kk,:),'r*',xx,vnum(kk,:),'b-')
    pause(0.0001)
end
dt=deltat;
t=(0:dt:final_t);
N=length(t);
plot(t,unum1(1:N),'r*',t,vnum1(1:N),'b-')
%axis([0,7,0,0.5]);
title(' solution at the column outlet')
legend('numerical solution c_1','numerical solution c_2')
xlabel('t (min)');
ylabel('c (g/l)');

% 使用 meshgrid 将 t 和 xx 转换为网格形式
[X, T] = meshgrid(xx, t);

figure;
subplot(1,2,1);
h1 = surf(X, T, unum(1:end-1,:));            % 用第四个参数默认就是 Z 本身来着色
set(h1, 'EdgeColor', 'none');     % 去掉网格线
shading interp                     % 平滑插值着色
colormap(jet)                      % 选用 jet 调色板
colorbar                           % 显示色条
title('numerical solution c_1');
xlabel('x');
ylabel('t');
zlabel('c_1');

subplot(1,2,2);
h2 = surf(X, T, vnum(1:end-1,:));           % 假设 unum2 是 c2 的数值解
set(h2, 'EdgeColor', 'none');
shading interp
colormap(parula)                  % 也可以试试不同调色板
colorbar
title('numerical solution c_2');
xlabel('x');
ylabel('t');
zlabel('c_2');

