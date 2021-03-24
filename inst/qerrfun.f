c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c       this is the end of the debugging code and the beginning
c       of the error function routines
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c 
c 
c 
        subroutine qerrfun(x7,erfun)
        implicit real *8 (a-h,o-z)
        save
        dimension coefs(18),coefs2(13),bufcoe(2)
        equivalence (coefs(18),bufcoe(1)),(coefs2(1),bufcoe(2))
c 
c        this subroutine evaluates the error function defined
c        by the formula
c 
c        erf(x)=2 /\pi * \int_0^x exp(-t**2) dt
c 
c       with extended precision accuracy.
c 
c 
c                      input parameters:
c 
c  x7 - the argument
c                      output parameters:
c 
c  erfun - error function of x
c 
        data coefs/
     1     0.46963291789115009827412922472047D+00,
     2     0.42582445804381043569204735290848D+00,
     3     -.49552626796204340403576830798593D-01,
     4     -.44929348876838274955800124201535D-02,
     5     0.12919410465849695349422476112098D-02,
     6     0.18363892921493962704169788672122D-04,
     7     -.22111147040995262915385554699197D-04,
     8     0.52333748523425713467369319779734D-06,
     9     0.27818478883353788538253099109576D-06,
     a     -.14115809274881311456031659199125D-07,
     1     -.27257129633056169998453886261514D-08,
     2     0.20634390487207062940639851585784D-09,
     3     0.21427399199678536792416443571655D-10,
     4     -.22299025553935820457965847853047D-11,
     5     -.13625007465069828057004944758513D-12,
     6     0.19514401092229309192942348364070D-13,
     7     0.68562716923170458929126778282939D-15,
     8     -.14450649286969997803041431643477D-15/
         data coefs2/
     9     -.24593530646054176842074301775172D-17,
     a     0.92959956122049262809396442021631D-18,
     1     0.28743670971555272132905892954155D-20,
     2     -.52871778042148820515979592731087D-20,
     3     0.44807385096647751700661175657819D-22,
     4     0.26922569397181593320864331505555D-22,
     5     -.50472040668074925808180499798587D-24,
     6     -.12386920454438148903803863927236D-24,
     7     0.35241362346072537624979164129445D-26,
     8     0.51841341153297714329403112468134D-27,
     9     -.19798231818357131118309229653530D-28,
     a     -.19862928018160851717691808512573D-29,
     1     0.64298976752690731987886063592519D-31/
c 
        x=abs(x7)
c 
c        if x .gt. 0.1 - evaluate erf via erfc
c 
        if(x .lt. 1.0d0) goto 1400
        call qerrfunc(x,erfunc)
        erfun=1-erfunc
        if(x7 .lt. 0) erfun=-erfun
        return
 1400 continue
c 
c        x .lt. 0.1. use chebychev series to evaluate erfun
c 
        t=x*2-1
        n=30
        call qerrexev(coefs,n,t,erfun)
        if(x7 .lt. 0) erfun=-erfun
        return
        end
c 
c 
c 
c 
c 
        subroutine qerrfunc(x7,erfunc)
        implicit real *8 (a-h,o-z)
        save
        dimension c02n1(10),c02n2(10),c02n3(10),c02n4(8),
     1     c24n1(10),c24n2(10),c24n3(12),
     2     c46n1(10),c46n2(10),c46n3(8),
     3     c68n1(10),c68n2(10),c68n3(6),
     4      c810n1(10),c810n2(15),
ccc     5     c02(30),c24(30),c46(30),c68(30),c810(30)
     5     c02(31),c24(30),c46(30),c68(30),c810(30)
        equivalence
     1   (c02(1),c02n1(1)),(c02(11),c02n2(1)),(c02(21),c02n3(1)),
     2   (c24(1),c24n1(1)),(c24(11),c24n2(1)),(c24(21),c24n3(1)),
     3   (c46(1),c46n1(1)),(c46(11),c46n2(1)),(c46(21),c46n3(1)),
     4   (c68(1),c68n1(1)),(c68(11),c68n2(1)),(c68(21),c68n3(1)),
     5   (c810(1),c810n1(1)),(c810(11),c810n2(1)),
     6   (c02(31),c02n4(1))
c 
c        this subroutine evaluates the second error function defined
c        by the formula
c 
c        erfc(x)=2/\pi*\int_x^{\infty}  exp(-t**2) dt = 1 - erf(x)
c 
c       with extended precision accuracy
c 
c 
c                      input parameters:
c 
c  x7 - the argument
c                      output parameters:
c 
c  erfunc - second error function of x
  
c 
c   coefficients for the interval [0, 2]
c 
        data c02n1/
     1     0.52136087049960619781849508312991D+00,
     2     -.34477074547252352449342439419592D+00,
     3     0.99759779268200963479718559480923D-01,
     4     -.26063665256478758984350198490020D-01,
     5     0.62684634943190106999675369238913D-02,
     6     -.14061223858870650278211581416723D-02,
     7     0.29698557246548996692018339283689D-03,
     8     -.59485553901067240043910010542263D-04,
     9     0.11362960761135521370839946482502D-04,
     a     -.20794078634147761077310996754027D-05/
c 
        data c02n2/
  
     1     0.36590821201669711438642556301882D-06,
     2     -.62107261363621591843172860865957D-07,
     3     0.10195316699270051335580482106724D-07,
     4     -.16223227802469610796062165768710D-08,
     5     0.25073406746662080536115811286194D-09,
     6     -.37703813593742712974168897019594D-10,
     7     0.55248768982320522072097020012181D-11,
     8     -.78999129790218407822939938977739D-12,
     9     0.11036213151407694552179330257083D-12,
     a     -.15079931632503800581334842519161D-13/
c 
        data c02n3/
     1     0.20174311859216235287859448030821D-14,
     2     -.26449521072902580164793348211251D-15,
     3     0.34011233662469910562850493803972D-16,
     4     -.42928773664251293839650309006466D-17,
     5     0.53223804776089466438685653147602D-18,
     6     -.64860644065140435856158116762429D-19,
     7     0.77739252371807837876370560162840D-20,
     8     -.91691812750698378690597135480207D-21,
     9     0.10648355585514911368237797170254D-21,
     a     -.12181814947819617722300168602698D-22/
c 
        data c02n4/
     1     0.13734755945084779867462110216846D-23,
     2     -.15268575573384197852242343081360D-24,
     3     0.16742427117522907125140060745447D-25,
     4     -.18114692950611147766479753260373D-26,
     5     0.19349079790342965313845605002607D-27,
     6     -.20546472152712031825029589284034D-28,
     7     0.23436671585903654146398399582971D-29,
     8     -.51064772345483477151570712662071D-30/
c 
c   coefficients for the interval [2, 4]
        data c24n1/
     1     0.18742956902622409306925144479872D+00,
     2     -.57946528836690967785413378444440D-01,
     3     0.85952473873965326996308911399058D-02,
     4     -.12284613640242593736724370562102D-02,
     5     0.16974366686180486953367670784488D-03,
     6     -.22738329337044562795469514814350D-04,
     7     0.29598443788048780996722287708018D-05,
     8     -.37513575581379711371503352440765D-06,
     9     0.46372828761738247083074639709199D-07,
     a     -.55994553925640879241202334507183D-08/
c 
        data c24n2/
     1     0.66131608966799159755426938686925D-09,
     2     -.76482715221998167746578867813974D-10,
     3     0.86709043258765937469529030691067D-11,
     4     -.96454562743821085144258450212591D-12,
     5     0.10536827959360625660700528989713D-12,
     6     -.11312630906760667269865441308807D-13,
     7     0.11945180669174576517515458382592D-14,
     8     -.12413107525804078010209734851108D-15,
     9     0.12702463231854429281216324752301D-16,
     a     -.12807258736639517232877482476199D-17/
c 
        data c24n3/
     1     0.12729426443963546096673617369866D-18,
     2     -.12478254325373568523046115484300D-19,
     3     0.12069382415080905790174026421391D-20,
     4     -.11523481357557426836576718890519D-21,
     5     0.10864746405132768883292074304440D-22,
     6     -.10119339333958407890345558606455D-23,
     7     0.93138978055373967457176232654886D-25,
     8     -.84741890625381197262746465981092D-26,
     9     0.76236749091460841177609634990828D-27,
     a     -.67781247299713562765608554360546D-28,
     1     0.58847289741147701455688797529953D-29,
     2     -.43217491793857309410573595035914D-30/
c 
c   coefficients for the interval [4, 6]
c 
        data c46n1/
     1     0.11277819743166245319088851431478D+00,
     2     -.21913497155272348854403231709269D-01,
     3     0.20915397845687494704395083143931D-02,
     4     -.19628624819122982453937138730095D-03,
     5     0.18126654238718201380392151756897D-04,
     6     -.16483628246572165758417591008001D-05,
     7     0.14769699327777976107481836575099D-06,
     8     -.13047490253442587780462366837498D-07,
     9     0.11369753471043668813639541094894D-08,
     a     -.97781802266404980741610121580943D-10/
c 
        data c46n2/
     1     0.83032123480629757872312378628952D-11,
     2     -.69646508425907134583559420164119D-12,
     3     0.57728321648271661335653788247210D-13,
     4     -.47301485583698741521430457595131D-14,
     5     0.38327149744119831383914884472204D-15,
     6     -.30720233149987745126274517435573D-16,
     7     0.24364677192448132179513790068467D-17,
     8     -.19126665070736228731222529219220D-18,
     9     0.14865432561241565654130238348592D-19,
     a     -.11441591030872881659849514692371D-20/
c 
        data c46n3/
     1     0.87230998571718318806405344232366D-22,
     2     -.65891651936999413334501938612018D-23,
     3     0.49324158413158222073023864676717D-24,
     4     -.36597300135688308474178347999623D-25,
     5     0.26920767745242577505300186270589D-26,
     6     -.19637445933108154996248853884989D-27,
     7     0.14216229691161037924544550733445D-28,
     8     -.10136635733468428578316773906482D-29/
c 
c   coefficients for the interval [6, 8]
c 
        data c68n1/
     1     0.80586726359943425276197738287254D-01,
     2     -.11340871814321887442002116443655D-01,
     3     0.79038930590881203560170174393107D-03,
     4     -.54574581010721115202842178559706D-04,
     5     0.37342298952738347988298869910226D-05,
     6     -.25326406802593954485737513767179D-06,
     7     0.17029548888897999275400134926829D-07,
     8     -.11354810726045605056112156438796D-08,
     9     0.75091260369026600390440269335786D-10,
     a     -.49262100700803286822366796153529D-11/
c 
        data c68n2/
     1     0.32064751085119214716490647595778D-12,
     2     -.20711328658873924858613392100474D-13,
     3     0.13277721023908436484100441115886D-14,
     4     -.84497101884086011030267767173465D-16,
     5     0.53385959541556925857276665429843D-17,
     6     -.33491977405720461950783951234504D-18,
     7     0.20866134081155284916762621651560D-19,
     8     -.12911813267013401127435041289982D-20,
     9     0.79365332146565188043514089609902D-22,
     a     -.48464641692687441074256035233402D-23/
c 
        data c68n3/
     1     0.29404919865774879091294502596371D-24,
     2     -.17728209305496042549085176793617D-25,
     3     0.10622256869982957270388264515729D-26,
     4     -.63281696042167086852717865443297D-28,
     5     0.37634510526886755284206278477043D-29,
     6     -.23276123507863826776390035921707D-30/
c 
c   coefficients for the interval [8,10]
c 
        data c810n1/
     1     0.62684290243453396610825314277179D-01,
     2     -.69014792436725362656233676014409D-02,
     3     0.37767451908518504573804627058453D-03,
     4     -.20547526563664360267957542649865D-04,
     5     0.11115026063587423010624134024677D-05,
     6     -.59787671676642149912916275307666D-07,
     7     0.31981785361812047436886673174224D-08,
     8     -.17014609183903509133469574199210D-09,
     9     0.90033981117773496809677415990439D-11,
     a     -.47390564596554209182714435555536D-12/
c 
        data c810n2/
     1     0.24814917515262421373897053355326D-13,
     2     -.12927149377630089157277771898625D-14,
     3     0.67002972077760909937376758563788D-16,
     4     -.34555574320072334655326037049359D-17,
     5     0.17733943812920036876032543058533D-18,
     6     -.90570174500582326678144704124431D-20,
     7     0.46034811426248772172003339213424D-21,
     8     -.23288236840799981234004091731405D-22,
     9     0.11726351904475800826720593142939D-23,
     a     -.58774876038377951007211157322197D-25,
     1     0.29325940462141413026313215565472D-26,
     2     -.14569212683159290298634440553278D-27,
     3     0.72274922213980504030355407245095D-29,
     4     -.37605098717833700992926894313892D-30,
     5     0.31974386737846878046756008801675D-31/
c 
        x=abs(x7)
c 
c       if x .gt. 10 - set the function to zero and exit
c 
        if(x .lt. 9.9999d0) goto 1100
        erfunc=0
        return
 1100 continue
c 
c        if x is on the interval [0,2] - act accordingly
c 
        if(x .gt. 2) goto 1200
        n=37
        t=x-1
        call qerrexev(c02,n,t,erfunc)
        erfunc=erfunc*exp(-x**2)
        if(x7 .lt. 0) erfunc=2-erfunc
        return
 1200 continue
c 
        if(x .gt. 4) goto 1400
        n=31
        t=x-3
        call qerrexev(c24,n,t,erfunc)
        erfunc=erfunc*exp(-x**2)
        if(x7 .lt. 0) erfunc=2-erfunc
        return
 1400 continue
c 
        if(x .gt. 6) goto 1600
        n=27
        t=x-5
        call qerrexev(c46,n,t,erfunc)
        erfunc=erfunc*exp(-x**2)
        if(x7 .lt. 0) erfunc=2-erfunc
        return
 1600 continue
c 
        if(x .gt. 8) goto 1800
        n=25
        t=x-7
        call qerrexev(c68,n,t,erfunc)
        erfunc=erfunc*exp(-x**2)
        if(x7 .lt. 0) erfunc=2-erfunc
        return
 1800 continue
c 
        if(x .gt. 10) goto 2000
        n=24
        t=x-9
        call qerrexev(c810,n,t,erfunc)
        erfunc=erfunc*exp(-x**2)
        if(x7 .lt. 0) erfunc=2-erfunc
        return
 2000 continue
  
  
  
        return
        end
c 
c 
c 
c 
c 
        subroutine qerrexev(c,n,x,fun)
        implicit real *8 (a-h,o-z)
        save
        dimension c(1)
c 
c        evaluate the chebychev expansion with coefficients
c        c at the point x
c 
cccc        call prinq('in qerrexev, x=*',x,1)
cccc        call prinq('in qerrexev, c=*',c,10)
        x2=x+x
        tkm1=1
        tk=x
        fun=c(1)*tkm1
        do 1400 k=1,n
        tkp1=x2*tk-tkm1
c 
        fun=fun+c(k+1)*tk
c 
        tkm1=tk
        tk=tkp1
 1400 continue
cccc        call prinq('and fun=*',fun,1)
        return
        end
