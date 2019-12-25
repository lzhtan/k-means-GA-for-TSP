function varargout =tsp_ga(u,re,m,n1,n2,n3,date1,date2,date3,dmat1,dmat2,dmat3,popSize,numIter,showProg,showResult);





% Sanity Checks
popSize = 4*ceil(popSize/4);
numIter = max(1,round(real(numIter(1))));
showProg = logical(showProg(1));
showResult = logical(showResult(1));
dims=2;
% Initialize the Population 初始化种群
pop1 = zeros(popSize,n1); %pop为500*200的矩阵
pop2 = zeros(popSize,n2); %pop为500*200的矩阵
pop3 = zeros(popSize,n3); %pop为500*200的矩阵
pop1(1,:) = (1:n1);%pop第一行填充“1、2、3、4.....499、500”
for k = 2:popSize
    pop1(k,:) = randperm(n1);%pop第其余行随机填充1——500之间的数字
end
pop2(1,:) = (1:n2);%pop第一行填充“1、2、3、4.....499、500”
for k = 2:popSize
    pop2(k,:) = randperm(n2);%pop第其余行随机填充1——500之间的数字
end
pop3(1,:) = (1:n3);%pop第一行填充“1、2、3、4.....499、500”
for k = 2:popSize
    pop3(k,:) = randperm(n3);%pop第其余行随机填充1——500之间的数字
end
% Run the GA
globalMin1 = Inf;  %globalMin为中间变量，记录全局最优路径长度，初始化为无穷大。
globalMin2 = Inf; 
globalMin3 = Inf; 
totalDist1 = zeros(1,popSize); %totalDist为中间变量，记录种群中500个个体各自的路径总长度。
globalMin2 = Inf; 
globalMin3 = Inf; 
distHistory1 = zeros(1,numIter);%distHistory为中间变量，记录迭代历史次数中globalMin。
distHistory2 = zeros(1,numIter);%distHistory为中间变量，记录迭代历史次数中globalMin。
distHistory3 = zeros(1,numIter);%distHistory为中间变量，记录迭代历史次数中globalMin。
tmpPop1 = zeros(4,n1);%tmpPop在变异过程中记录挑选的4组父母的具体路径。
newPop1 = zeros(popSize,n1);%newPop在变异过程中记录变异后子代的具体路径，为500*200矩阵。
tmpPop2 = zeros(4,n2);%tmpPop在变异过程中记录挑选的4组父母的具体路径。
newPop2 = zeros(popSize,n2);%newPop在变异过程中记录变异后子代的具体路径，为500*200矩阵。
tmpPop3 = zeros(4,n3);%tmpPop在变异过程中记录挑选的4组父母的具体路径。
newPop3 = zeros(popSize,n3);%newPop在变异过程中记录变异后子代的具体路径，为500*200矩阵。
if showProg
    pfig = figure('Name','TSP_GA | Current Best Solution','Numbertitle','off');
end
for iter = 1:numIter %迭代次数
    % Evaluate Each Population Member (Calculate Total Distance)
    for p = 1:popSize
        d = dmat1(pop1(p,n1),pop1(p,1)); % Closed Path 选择算子
        for k = 2:n1
            d = d + dmat1(pop1(p,k-1),pop1(p,k));  %交叉算子、变异算子
        end
        totalDist1(p) = d;  % 计算适应度 
    end
        for p = 1:popSize
        d = dmat3(pop3(p,n3),pop3(p,1)); % Closed Path 选择算子
        for k = 2:n3
            d = d + dmat3(pop3(p,k-1),pop3(p,k));  %交叉算子、变异算子
        end
        totalDist3(p) = d;  % 计算适应度 
        end
        for p = 1:popSize
        d = dmat2(pop2(p,n2),pop2(p,1)); % Closed Path 选择算子
        for k = 2:n2
            d = d + dmat2(pop2(p,k-1),pop2(p,k));  %交叉算子、变异算子
        end
        totalDist2(p) = d;  % 计算适应度 
    end
    
    %找到最小和最大适应度的染色体及它们在种群中的位置
    % Find the Best Route in the Population
    [minDist1,index] = min(totalDist1);
    distHistory1(iter) = minDist1;
     % 代替上一次进化中最好的染色体
    if minDist1 < globalMin1
        globalMin1 = minDist1;
        optRoute1 = pop1(index,:);
        if showProg
            % Plot the Best Route
            figure(pfig);
            rte1 = optRoute1([1:n1 1]);
            if dims > 2, 
                plot3(date1(rte1,1),date1(rte1,2),date1(rte1,3),'r.-');
            else
                plot(date1(rte1,1),date1(rte1,2),'r.-');
            end
            title(sprintf('Total Distance = %1.4f, Iteration = %d',minDist1,iter));
        end
    end
    [minDist2,index] = min(totalDist2);
    distHistory2(iter) = minDist2;
     % 代替上一次进化中最好的染色体
    if minDist2 < globalMin2
        globalMin2 = minDist2;
        optRoute2 = pop2(index,:);
        if showProg
            % Plot the Best Route
            figure(pfig);
            rte2 = optRoute2([1:n2 1]);
            if dims > 2, 
                plot3(date2(rte2,1),date2(rte2,2),date2(rte2,3),'r.-');
            else
                plot(date2(rte2,1),date2(rte2,2),'r.-');
            end
            title(sprintf('Total Distance = %1.4f, Iteration = %d',minDist2,iter));
        end
    end
        [minDist3,index] = min(totalDist3);
    distHistory3(iter) = minDist3;
     % 代替上一次进化中最好的染色体
    if minDist3 < globalMin3
        globalMin3 = minDist3;
        optRoute3 = pop3(index,:);
   %     if showProg
    %        % Plot the Best Route
   %         figure(pfig);
   %         rte3 = optRoute3([1:n3 1]);
   %         if dims > 2, 
   %             plot3(date3(rte3,1),date3(rte3,2),date3(rte3,3),'r.-');
   %         else
  %              plot(date3(rte3,1),date3(rte3,2),'r.-');
    %        end
   %         title(sprintf('Total Distance = %1.4f, Iteration = %d',minDist3,iter));
  %  end
        
    end
    
    % Genetic Algorithm Operators
    randomOrder = randperm(popSize);%randperm是matlab函数，功能是随机打乱一个数字序列。　　函数功能：随机打乱一个数字序列。　　语法格式：　　y = randperm(n)　　y是把1到n这些数随机打乱得到的一个数字序列。
    %randomOrder为打乱顺序的种群总数。级从0-500打乱顺序。
    for p = 4:4:popSize    %每间隔四个跳一步，从0-500跳完。
        rtes1 = pop1(randomOrder(p-3:p),:);%rtes记录取出的四组父母。为4*200矩阵。
        dists1 = totalDist1(randomOrder(p-3:p));  %dists记录四组父母的路径总距离，为4个值。
        [ignore,idx] = min(dists1); %#ok %[ignore,idx]储存min(dists)
        bestOf4Route1 = rtes1(idx,:);%bestOf4Route储存当前四组父母中最优具体路径
        routeInsertionPoints1 = sort(ceil(n1*rand(1,2)));%rand(1,2)随机出两个0-1的数，例如0.0754、0.2734；ceil向上取整；sort将取出的城市数从小到大排序，然后给routeInsertionPoints。
        I = routeInsertionPoints1(1);
        J = routeInsertionPoints1(2);%将routeInsertionPoints中的两个数分别给I和J
        for k = 1:4 % Mutate the Best to get Three New Routes  交叉产生新的个体  四组父母中取最优个体进行不同操作并将其结果覆盖原较差三组。
            tmpPop1(k,:) = bestOf4Route1;
            switch k
                case 2 % Flip I-J列交换
                    tmpPop1(k,I:J) = tmpPop1(k,J:-1:I);
                case 3 % Swap 只I、J两列交换
                    tmpPop1(k,[I J]) = tmpPop1(k,[J I]);
                case 4 % Slide
                    tmpPop1(k,I:J) = tmpPop1(k,[I+1:J I]);
                otherwise % Do Nothing
            end
        end
        newPop1(p-3:p,:) = tmpPop1;%当前四组放入子代
    end
    pop1 = newPop1;


% Genetic Algorithm Operators
    randomOrder = randperm(popSize);%randperm是matlab函数，功能是随机打乱一个数字序列。　　函数功能：随机打乱一个数字序列。　　语法格式：　　y = randperm(n)　　y是把1到n这些数随机打乱得到的一个数字序列。
    %randomOrder为打乱顺序的种群总数。级从0-500打乱顺序。
    for p = 4:4:popSize    %每间隔四个跳一步，从0-500跳完。
        rtes2 = pop2(randomOrder(p-3:p),:);%rtes记录取出的四组父母。为4*200矩阵。
        dists2 = totalDist2(randomOrder(p-3:p));  %dists记录四组父母的路径总距离，为4个值。
        [ignore,idx] = min(dists2); %#ok %[ignore,idx]储存min(dists)
        bestOf4Route2 = rtes2(idx,:);%bestOf4Route储存当前四组父母中最优具体路径
        routeInsertionPoints2 = sort(ceil(n2*rand(1,2)));%rand(1,2)随机出两个0-1的数，例如0.0754、0.2734；ceil向上取整；sort将取出的城市数从小到大排序，然后给routeInsertionPoints。
        I = routeInsertionPoints2(1);
        J = routeInsertionPoints2(2);%将routeInsertionPoints中的两个数分别给I和J
        for k = 1:4 % Mutate the Best to get Three New Routes  交叉产生新的个体  四组父母中取最优个体进行不同操作并将其结果覆盖原较差三组。
            tmpPop2(k,:) = bestOf4Route2;
            switch k
                case 2 % Flip I-J列交换
                    tmpPop2(k,I:J) = tmpPop2(k,J:-1:I);
                case 3 % Swap 只I、J两列交换
                    tmpPop2(k,[I J]) = tmpPop2(k,[J I]);
                case 4 % Slide
                    tmpPop2(k,I:J) = tmpPop2(k,[I+1:J I]);
                otherwise % Do Nothing
            end
        end
        newPop2(p-3:p,:) = tmpPop2;%当前四组放入子代
    end
    pop2 = newPop2;


% Genetic Algorithm Operators
    randomOrder = randperm(popSize);%randperm是matlab函数，功能是随机打乱一个数字序列。　　函数功能：随机打乱一个数字序列。　　语法格式：　　y = randperm(n)　　y是把1到n这些数随机打乱得到的一个数字序列。
    %randomOrder为打乱顺序的种群总数。级从0-500打乱顺序。
    for p = 4:4:popSize    %每间隔四个跳一步，从0-500跳完。
        rtes3 = pop3(randomOrder(p-3:p),:);%rtes记录取出的四组父母。为4*200矩阵。
        dists3 = totalDist3(randomOrder(p-3:p));  %dists记录四组父母的路径总距离，为4个值。
        [ignore,idx] = min(dists3); %#ok %[ignore,idx]储存min(dists)
        bestOf4Route3 = rtes3(idx,:);%bestOf4Route储存当前四组父母中最优具体路径
        routeInsertionPoints3 = sort(ceil(n3*rand(1,2)));%rand(1,2)随机出两个0-1的数，例如0.0754、0.2734；ceil向上取整；sort将取出的城市数从小到大排序，然后给routeInsertionPoints。
        I = routeInsertionPoints3(1);
        J = routeInsertionPoints3(2);%将routeInsertionPoints中的两个数分别给I和J
        for k = 1:4 % Mutate the Best to get Three New Routes  交叉产生新的个体  四组父母中取最优个体进行不同操作并将其结果覆盖原较差三组。
            tmpPop3(k,:) = bestOf4Route3;
            switch k
                case 2 % Flip I-J列交换
                    tmpPop3(k,I:J) = tmpPop3(k,J:-1:I);
                case 3 % Swap 只I、J两列交换
                    tmpPop3(k,[I J]) = tmpPop3(k,[J I]);
                case 4 % Slide
                    tmpPop3(k,I:J) = tmpPop3(k,[I+1:J I]);
                otherwise 
                     % Do Nothing
            end
        end
        newPop3(p-3:p,:) = tmpPop3;%当前四组放入子代
    end
    pop3 = newPop3;
end



%if (u(1,1)<u(2,1)&&u(2,1)<u(3,1))
d = zeros(8); 
bijiao=Inf;
for k = 2:(n1-1)
    for kk = 2:(n2-1)
    d(1)=sqrt((date1(optRoute1(k),1)-date2(optRoute2(kk-1),1))^2+(date1(optRoute1(k),2)-date2(optRoute2(kk-1),2))^2)+sqrt((date1(optRoute1(k-1),1)-date2(optRoute2(kk),1))^2+(date1(optRoute1(k-1),2)-date2(optRoute2(kk),2))^2)-sqrt((date1(optRoute1(k),1)-date1(optRoute1(k-1),1))^2+(date1(optRoute1(k),2)-date1(optRoute1(k-1),2))^2)-sqrt((date2(optRoute2(kk),1)-date2(optRoute2(kk-1),1))^2+(date2(optRoute2(kk),1)-date2(optRoute2(kk-1),1))^2);
    if d(1)< bijiao
       bijiao=d(1);
       jiaochadian1 = k;jiaochadian2 = k-1,jiaochadian3=kk-1,jiaochadian4=kk;
   end
    d(2)=sqrt((date1(optRoute1(k-1),1)-date2(optRoute2(kk-1),1))^2+(date1(optRoute1(k-1),2)-date2(optRoute2(kk-1),2))^2)+sqrt((date1(optRoute1(k),1)-date2(optRoute2(kk),1))^2+(date1(optRoute1(k),2)-date2(optRoute2(kk),2))^2)-sqrt((date1(optRoute1(k),1)-date1(optRoute1(k-1),1))^2+(date1(optRoute1(k),2)-date1(optRoute1(k-1),2))^2)-sqrt((date2(optRoute2(kk),1)-date2(optRoute2(kk-1),1))^2+(date2(optRoute2(kk),1)-date2(optRoute2(kk-1),1))^2);
   if d(2)< bijiao
       bijiao=d(2);
       jiaochadian1=k-1,jiaochadian2=k,jiaochadian3=kk-1,jiaochadian4=kk;
   end
   d(3)=sqrt((date1(optRoute1(k),1)-date2(optRoute2(kk-1),1))^2+(date1(optRoute1(k),2)-date2(optRoute2(kk-1),2))^2)+sqrt((date1(optRoute1(k+1),1)-date2(optRoute2(kk),1))^2+(date1(optRoute1(k+1),2)-date2(optRoute2(kk),2))^2)-sqrt((date1(optRoute1(k),1)-date1(optRoute1(k+1),1))^2+(date1(optRoute1(k),2)-date1(optRoute1(k+1),2))^2)-sqrt((date2(optRoute2(kk),1)-date2(optRoute2(kk-1),1))^2+(date2(optRoute2(kk),1)-date2(optRoute2(kk-1),1))^2);
     if d(3)< bijiao
        bijiao=d(3);
       jiaochadian1=k,jiaochadian2=k+1,jiaochadian3=kk-1,jiaochadian4=kk;
     end
     d(4)=sqrt((date1(optRoute1(k+1),1)-date2(optRoute2(kk-1),1))^2+(date1(optRoute1(k+1),2)-date2(optRoute2(kk-1),2))^2)+sqrt((date1(optRoute1(k),1)-date2(optRoute2(kk),1))^2+(date1(optRoute1(k),2)-date2(optRoute2(kk),2))^2)-sqrt((date1(optRoute1(k),1)-date1(optRoute1(k+1),1))^2+(date1(optRoute1(k),2)-date1(optRoute1(k+1),2))^2)-sqrt((date2(optRoute2(kk),1)-date2(optRoute2(kk-1),1))^2+(date2(optRoute2(kk),1)-date2(optRoute2(kk-1),1))^2);
     if d(4)< bijiao
       bijiao=d(4);
       jiaochadian1=k+1,jiaochadian2=k,jiaochadian3=kk-1,jiaochadian4=kk;
     end
   d(5)=sqrt((date1(optRoute1(k),1)-date2(optRoute2(kk),1))^2+(date1(optRoute1(k),2)-date2(optRoute2(kk),2))^2)+sqrt((date1(optRoute1(k-1),1)-date2(optRoute2(kk+1),1))^2+(date1(optRoute1(k-1),2)-date2(optRoute2(kk+1),2))^2)-sqrt((date1(optRoute1(k),1)-date1(optRoute1(k-1),1))^2+(date1(optRoute1(k),2)-date1(optRoute1(k-1),2))^2)-sqrt((date2(optRoute2(kk),1)-date2(optRoute2(kk+1),1))^2+(date2(optRoute2(kk),1)-date2(optRoute2(kk+1),1))^2);
   if d(5)< bijiao
        bijiao=d(5);
       jiaochadian1=k,jiaochadian2=k-1,jiaochadian3=kk,jiaochadian4=kk+1;
   end
   d(6)=sqrt((date1(optRoute1(k),1)-date2(optRoute2(kk+1),1))^2+(date1(optRoute1(k),2)-date2(optRoute2(kk+1),2))^2)+sqrt((date1(optRoute1(k-1),1)-date2(optRoute2(kk),1))^2+(date1(optRoute1(k-1),2)-date2(optRoute2(kk),2))^2)-sqrt((date1(optRoute1(k),1)-date1(optRoute1(k-1),1))^2+(date1(optRoute1(k),2)-date1(optRoute1(k-1),2))^2)-sqrt((date2(optRoute2(kk),1)-date2(optRoute2(kk+1),1))^2+(date2(optRoute2(kk),1)-date2(optRoute2(kk+1),1))^2);
      if d(6)< bijiao
        bijiao=d(6);
       jiaochadian1=k,jiaochadian2=k-1,jiaochadian3=kk+1,jiaochadian4=kk;
    end
   d(7)=sqrt((date1(optRoute1(k),1)-date2(optRoute2(kk+1),1))^2+(date1(optRoute1(k),2)-date2(optRoute2(kk+1),2))^2)+sqrt((date1(optRoute1(k+1),1)-date2(optRoute2(kk),1))^2+(date1(optRoute1(k+1),2)-date2(optRoute2(kk),2))^2)-sqrt((date1(optRoute1(k),1)-date1(optRoute1(k+1),1))^2+(date1(optRoute1(k),2)-date1(optRoute1(k+1),2))^2)-sqrt((date2(optRoute2(kk),1)-date2(optRoute2(kk+1),1))^2+(date2(optRoute2(kk),1)-date2(optRoute2(kk+1),1))^2);
     if d(7)< bijiao
       bijiao=d(7);
       jiaochadian1=k,jiaochadian2=k+1,jiaochadian3=kk+1,jiaochadian4=kk;
     end
     d(8)=sqrt((date1(optRoute1(k),1)-date2(optRoute2(kk),1))^2+(date1(optRoute1(k),2)-date2(optRoute2(kk),2))^2)+sqrt((date1(optRoute1(k+1),1)-date2(optRoute2(kk+1),1))^2+(date1(optRoute1(k+1),2)-date2(optRoute2(kk+1),2))^2)-sqrt((date1(optRoute1(k),1)-date1(optRoute1(k+1),1))^2+(date1(optRoute1(k),2)-date1(optRoute1(k+1),2))^2)-sqrt((date2(optRoute2(kk),1)-date2(optRoute2(kk+1),1))^2+(date2(optRoute2(kk),1)-date2(optRoute2(kk+1),1))^2);
    if d(8)< bijiao
          bijiao=d(8);
       jiaochadian1=k,jiaochadian2=k+1,jiaochadian3=kk,jiaochadian4=kk+1;
   end
    end
end
jiaochadian1
jiaochadian2
jiaochadian3
jiaochadian4

if jiaochadian1<jiaochadian2
for k = 1:jiaochadian1
    optRoutee(k,:)=date1(optRoute1(k),:);
end
 if jiaochadian3 < jiaochadian4%没问题
    for kk = jiaochadian3:-1:1
    optRoutee((jiaochadian1+jiaochadian3-kk+1),:)=date2(optRoute2(kk),:);
    end
    for  kkk = n2:-1:jiaochadian4
    optRoutee((jiaochadian1+jiaochadian3+n2-kkk+1),:)=date2(optRoute2(kkk),:);
    end
    for  kkkkk = jiaochadian2 : n1
    optRoutee((n2+kkkkk),:)=date1(optRoute1(kkkkk),:);
    end
 else     %没问题
    for k = jiaochadian3 : n2
    optRoutee((jiaochadian1+k-jiaochadian3+1),:)=date2(optRoute2(k),:);
    end
     for k = 1 : jiaochadian4
    optRoutee((jiaochadian1+n2-jiaochadian4+k),:)=date2(optRoute2(k),:);
     end
    for k = jiaochadian2 : n1
    optRoutee((n2+k),:)=date1(optRoute1(k),:);
    end
 end
end


if jiaochadian1>jiaochadian2
for k = 1:jiaochadian2
    optRoutee(k,:)=date1(optRoute1(k),:);
end
 if jiaochadian3 < jiaochadian4
    for k = jiaochadian4 : n2
    optRoutee((jiaochadian2+k-jiaochadian4+1),:)=date2(optRoute2(k),:);
    end
     for  k = 1 : jiaochadian3
    optRoutee((jiaochadian2+n2+k-jiaochadian2),:)=date2(optRoute2(k),:);
     end
     for  k = jiaochadian1 :-1: 1
    optRoutee((jiaochadian1+n2+jiaochadian1-k+1),:)=date1(optRoute1(k),:);
    end
  else
    for k = jiaochadian4:-1 : 1
        optRoutee((jiaochadian2+jiaochadian4-k+1),:)=date2(optRoute2(k),:);
    end
     for k = n2 : -1 :jiaochadian3
        optRoutee((jiaochadian2+jiaochadian4+n2-k+1),:)=date2(optRoute2(k),:);
     end
    for k = jiaochadian1 : n1
    optRoutee((n2+k),:)=date1(optRoute1(k),:);
    end
 end
end
%end




%if (u(1,1)<u(2,1)&&u(2,1)<u(3,1))
d = zeros(8); 
bijiao=Inf;
for k = 2:(n1+n2-1)
    for kk = 2:(n3-1)
    d(1)=sqrt((optRoutee((k),1)-date3(optRoute3(kk-1),1))^2+(optRoutee((k),2)-date3(optRoute3(kk-1),2))^2)+sqrt((optRoutee((k-1),1)-date3(optRoute3(kk),1))^2+(optRoutee((k-1),2)-date3(optRoute3(kk),2))^2)-sqrt((optRoutee((k),1)-optRoutee((k-1),1))^2+(optRoutee(k,2)-optRoutee((k-1),2))^2)-sqrt((date3(optRoute3(kk),1)-date3(optRoute3(kk-1),1))^2+(date3(optRoute3(kk),1)-date3(optRoute3(kk-1),1))^2);
    if d(1)< bijiao
       bijiao=d(1);
       jiaochadian1 = k;jiaochadian2 = k-1,jiaochadian3=kk-1,jiaochadian4=kk;
   end
    d(2)=sqrt((optRoutee((k-1),1)-date3(optRoute3(kk-1),1))^2+(optRoutee((k-1),2)-date3(optRoute3(kk-1),2))^2)+sqrt((optRoutee((k),1)-date3(optRoute3(kk),1))^2+(optRoutee((k),2)-date3(optRoute3(kk),2))^2)-sqrt((optRoutee((k),1)-optRoutee((k-1),1))^2+(optRoutee(k,2)-optRoutee((k-1),2))^2)-sqrt((date3(optRoute3(kk),1)-date3(optRoute3(kk-1),1))^2+(date3(optRoute3(kk),1)-date3(optRoute3(kk-1),1))^2);
    if d(2)< bijiao
       bijiao=d(2);
       jiaochadian1=k-1,jiaochadian2=k,jiaochadian3=kk-1,jiaochadian4=kk;
    end
    d(3)=sqrt((optRoutee((k),1)-date3(optRoute3(kk-1),1))^2+(optRoutee((k),2)-date3(optRoute3(kk-1),2))^2)+sqrt((optRoutee((k+1),1)-date3(optRoute3(kk),1))^2+(optRoutee((k+1),2)-date3(optRoute3(kk),2))^2)-sqrt((optRoutee((k),1)-optRoutee((k+1),1))^2+(optRoutee(k,2)-optRoutee((k+1),2))^2)-sqrt((date3(optRoute3(kk),1)-date3(optRoute3(kk-1),1))^2+(date3(optRoute3(kk),1)-date3(optRoute3(kk-1),1))^2);
     if d(3)< bijiao
        bijiao=d(3);
       jiaochadian1=k,jiaochadian2=k+1,jiaochadian3=kk-1,jiaochadian4=kk;
     end
    d(4)=sqrt((optRoutee((k+1),1)-date3(optRoute3(kk-1),1))^2+(optRoutee((k+1),2)-date3(optRoute3(kk-1),2))^2)+sqrt((optRoutee((k),1)-date3(optRoute3(kk),1))^2+(optRoutee((k),2)-date3(optRoute3(kk),2))^2)-sqrt((optRoutee((k),1)-optRoutee((k+1),1))^2+(optRoutee(k,2)-optRoutee((k+1),2))^2)-sqrt((date3(optRoute3(kk),1)-date3(optRoute3(kk-1),1))^2+(date3(optRoute3(kk),1)-date3(optRoute3(kk-1),1))^2);
     if d(4)< bijiao
       bijiao=d(4);
       jiaochadian1=k+1,jiaochadian2=k,jiaochadian3=kk-1,jiaochadian4=kk;
     end
    d(5)=sqrt((optRoutee((k),1)-date3(optRoute3(kk),1))^2+(optRoutee((k),2)-date3(optRoute3(kk),2))^2)+sqrt((optRoutee((k-1),1)-date3(optRoute3(kk+1),1))^2+(optRoutee((k-1),2)-date3(optRoute3(kk+1),2))^2)-sqrt((optRoutee((k),1)-optRoutee((k-1),1))^2+(optRoutee(k,2)-optRoutee((k-1),2))^2)-sqrt((date3(optRoute3(kk),1)-date3(optRoute3(kk+1),1))^2+(date3(optRoute3(kk),1)-date3(optRoute3(kk+1),1))^2);
     if d(5)< bijiao
        bijiao=d(5);
       jiaochadian1=k,jiaochadian2=k-1,jiaochadian3=kk,jiaochadian4=kk+1;
     end
    d(6)=sqrt((optRoutee((k),1)-date3(optRoute3(kk+1),1))^2+(optRoutee((k),2)-date3(optRoute3(kk+1),2))^2)+sqrt((optRoutee((k-1),1)-date3(optRoute3(kk),1))^2+(optRoutee((k-1),2)-date3(optRoute3(kk),2))^2)-sqrt((optRoutee((k),1)-optRoutee((k-1),1))^2+(optRoutee(k,2)-optRoutee((k-1),2))^2)-sqrt((date3(optRoute3(kk),1)-date3(optRoute3(kk+1),1))^2+(date3(optRoute3(kk),1)-date3(optRoute3(kk+1),1))^2);
     if d(6)< bijiao
        bijiao=d(6);
       jiaochadian1=k,jiaochadian2=k-1,jiaochadian3=kk+1,jiaochadian4=kk;
     end
    d(7)=sqrt((optRoutee((k),1)-date3(optRoute3(kk+1),1))^2+(optRoutee((k),2)-date3(optRoute3(kk+1),2))^2)+sqrt((optRoutee((k+1),1)-date3(optRoute3(kk),1))^2+(optRoutee((k+1),2)-date3(optRoute3(kk),2))^2)-sqrt((optRoutee((k),1)-optRoutee((k+1),1))^2+(optRoutee(k,2)-optRoutee((k+1),2))^2)-sqrt((date3(optRoute3(kk),1)-date3(optRoute3(kk+1),1))^2+(date3(optRoute3(kk),1)-date3(optRoute3(kk+1),1))^2);
   
      if d(7)< bijiao
       bijiao=d(7);
       jiaochadian1=k,jiaochadian2=k+1,jiaochadian3=kk+1,jiaochadian4=kk;
      end
    d(8)=sqrt((optRoutee((k),1)-date3(optRoute3(kk),1))^2+(optRoutee((k),2)-date3(optRoute3(kk),2))^2)+sqrt((optRoutee((k+1),1)-date3(optRoute3(kk+1),1))^2+(optRoutee((k+1),2)-date3(optRoute3(kk+1),2))^2)-sqrt((optRoutee((k),1)-optRoutee((k+1),1))^2+(optRoutee(k,2)-optRoutee((k+1),2))^2)-sqrt((date3(optRoute3(kk),1)-date3(optRoute3(kk+1),1))^2+(date3(optRoute3(kk),1)-date3(optRoute3(kk+1),1))^2);
     if d(8)< bijiao
          bijiao=d(8);
       jiaochadian1=k,jiaochadian2=k+1,jiaochadian3=kk,jiaochadian4=kk+1;
   end
    end
end


if jiaochadian1<jiaochadian2
for k = 1:jiaochadian1
    optRouteee(k,:)=optRoutee(k,:);
end
 if jiaochadian3 < jiaochadian4%没问题
    for kk = jiaochadian3:-1:1
    optRouteee((jiaochadian1+jiaochadian3-kk+1),:)=date3(optRoute3(kk),:);
    end
    for  kkk = n3:-1:jiaochadian4
    optRouteee((jiaochadian1+jiaochadian3+n3-kkk+1),:)=date3(optRoute3(kkk),:);
    end
    for  kkkkk = jiaochadian2 : (n1+n2)
    optRouteee((n3+kkkkk),:)=optRoutee(kkkkk,:);
    end
 else     %没问题
    for k = jiaochadian3 : n3
    optRouteee((jiaochadian1+k-jiaochadian3+1),:)=date3(optRoute3(k),:);
    end
     for k = 1 : jiaochadian4
    optRouteee((jiaochadian1+n3-jiaochadian4+k),:)=date3(optRoute3(k),:);
     end
    for k = jiaochadian2 : n1+n2
    optRouteee((n3+k),:)=optRoutee(k,:);
    end
 end
end


if jiaochadian1>jiaochadian2
for k = 1:jiaochadian2
    optRouteee(k,:)=optRoutee(k,:);
end
 if jiaochadian3 < jiaochadian4
    for k = jiaochadian4 : n3
    optRouteee((jiaochadian2+k-jiaochadian4+1),:)=date3(optRoute3(k),:);
    end
     for  k = 1 : jiaochadian3
    optRouteee((jiaochadian2+n3+k-jiaochadian3-jiaochadian4+1),:)=date3(optRoute3(k),:);
     end
     for  k = jiaochadian1 :(n1+n2)
    optRouteee((n3+k),:)=date1(optRoute1(k),:);
    end
  else
    for k = jiaochadian4:-1 : 1
        optRouteee((jiaochadian2+jiaochadian4-k+1),:)=date3(optRoute3(k),:);
    end
     for k = n3 : -1 :jiaochadian3
        optRouteee((jiaochadian2+jiaochadian4+n3-k+1),:)=date3(optRoute3(k),:);
     end
    for k = jiaochadian1 : n1+n2
    optRouteee((n3+k),:)=optRouteee(k,:);
    end
 end
end

d = sqrt((optRouteee(1,1)-(optRouteee(m,1)))^2+(optRouteee(1,2)-(optRouteee(m,2)))^2);
        for k = 2:m
            d = d + sqrt((optRouteee(k-1,1)-(optRouteee(k,1)))^2+(optRouteee(k-1,2)-(optRouteee(k,2)))^2);  %交叉算子、变异算子
        end

          figure(pfig);
          rte = 1:m;
          rte(m+1)=1;
          plot(optRouteee(rte,1),optRouteee(rte,2),'r.-');
          title(sprintf('Total Distance = %1.4f, Iteration = %d',d,iter));
   
% 画图
if showResult
    % Plots the GA Results
    figure('Name','TSP_GA | Results','Numbertitle','off');
    % subplot(2,2,1);
    figure(1),
    pclr = ~get(0,'DefaultAxesColor');
    if dims > 2, 
        plot3(date1(:,1),date1(:,2),date1(:,3),'.','Color',pclr);
    else
        plot(date1(:,1),date1(:,2),'.','Color',pclr)
        hold on
        plot(date3(:,1),date3(:,2),'.','Color',pclr)
        hold on
        plot(date3(:,1),date3(:,2),'.','Color',pclr)
    end
    title('城市位置');
    grid on
%     subplot(2,2,2);
    figure(3),
    rte1 = optRoute1([1:n1 1]);
    rte2 = optRoute2([1:n2 1]);
    rte3 = optRoute3([1:n3 1]);
    if dims > 2, plot3(date1(rte,1),date1(rte,2),date1(rte,3),'r.-');
    else
        plot(date1(rte1,1),date1(rte1,2),'r.-')
        hold on
        plot(date2(rte2,1),date2(rte2,2),'g.-')
        hold on
        plot(date3(rte3,1),date3(rte3,2),'b.-')
        hold on
    for i=1:m 
    if re(i,3)==1   
         plot(re(i,1),re(i,2),'ro'); 
    elseif re(i,3)==2
         plot(re(i,1),re(i,2),'go'); 
    else 
         plot(re(i,1),re(i,2),'bo'); 
    end
    end
    end
    title(sprintf('最短距离 = %1.4f',minDist1,minDist2,minDist3));
    grid on
%     subplot(2,2,4);
    for k=1:numIter
        distHistory(k)=distHistory1(k)+distHistory2(k)+distHistory3(k);
    end
    figure(4),
    plot(distHistory,'b','LineWidth',2);
    hold on
    title('Best fitness curve');
    grid on
    set(gca,'XLim',[0 numIter+1],'YLim',[0 1.1*max([1 distHistory])]);
end

% Return Outputs
if nargout
    varargout{1} = d;
    
end
