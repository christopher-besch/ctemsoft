        subroutine  eispack(nm,n,ar,ai,wr,wi,vr,vi,ierr)
c
c   eigenvalues and eigenvectors of a complex general matrix
c   calling EISPACK subroutine (B. T. Smith, Lecture notes in computer
c   science, Vol 6, Springer-Verlag 1976)
c   implemented by J.M. Zuo sept. 1987
c
        integer nm,n,is1,is2,iv1(n)
        real ar(nm,n),ai(nm,n),wr(n),wi(n),vr(nm,n),vi(nm,n),fv1(n)
        real fv2(n),fv3(n)
c
        call cbal(nm,n,ar,ai,is1,is2,fv1)
        call corth(nm,n,is1,is2,ar,ai,fv2,fv3)
        call comqr2(nm,n,is1,is2,fv2,fv3,ar,ai,wr,wi,vr,vi,ierr)
        if (ierr.ne.0) return
        call cbabk2(nm,n,is1,is2,fv1,n,vr,vi)
c
c   Normalize the eigenvectors
c
        do 160 i=1,n
          sum=0.0
          do 100 j=1,n
                term=vr(j,i)**2+vi(j,i)**2
                sum=sum+term
 100       continue
          sum=sqrt(sum)
          do 120 j=1,n
                vr(j,i)=vr(j,i)/sum
                vi(j,i)=vi(j,i)/sum
 120       continue
 160     continue
        return
        end


        subroutine cbal(nm,n,ar,ai,low,igh,scale)
c
        integer i,j,k,l,m,n,jj,nm,igh,low,iexc
        real ar(nm,n),ai(nm,n),scale(n)
        real c,f,g,r,s,b2,radix
        real abs
        logical noconv
        radix=2.0
        b2=radix*radix
        k=1
        l=n
        go to 100
20      scale(m)=j
        if (j.eq.m) goto 50
        do 30 i=1,l
          f=ar(i,j)
          ar(i,j)=ar(i,m)
          ar(i,m)=f
          f=ai(i,j)
          ai(i,j)=ai(i,m)
          ai(i,m)=f
30      continue
        do 40 i=k,n
          f=ar(j,i)
          ar(j,i)=ar(m,i)
          ar(m,i)=f
          f=ai(j,i)
          ai(j,i)=ai(m,i)
          ai(m,i)=f
40      continue
50      goto (80,130) iexc
80      if(l.eq.1) goto 280
        l=l-1
100     do 120 jj=1,l
          j=l+1-jj
          do 110 i=1,l
            if(i.eq.j) goto 110
            if(ar(j,i).ne.0.0.or.ai(j,i).ne.0.0) goto 120
110       continue
          m=l
          iexc=1
          go to 20
120     continue
        goto 140
130     k=k+1
140     do 170 j=k,l
          do 150 i=k,l
            if(i.eq.j) goto 150
            if(ar(i,j).ne.0.0.or.ai(i,j).ne.0.0) goto 170
150     continue
        m=k
        iexc=2
        goto 20
170     continue
        do 180 i=k,l
180       scale(i)=1.0
190     noconv=.false.
        do 270 i=k,l
          c=0.0
          r=0.0
          do 200 j=k,l
            if(j.eq.i) goto 200
            c=c+abs(ar(j,i))+abs(ai(j,i))
            r=r+abs(ar(i,j))+abs(ai(i,j))
200     continue
        g=r/radix
        f=1.0
        s=c+r
210     if(c.ge.g) goto 220
        f=f*radix
        c=c*b2
        goto 210
220     g=r*radix
230     if(c.lt.g) goto 240
        f=f/radix
        c=c/b2
        goto 230
240     if((c+r)/f.ge.0.95*s) goto 270
        g=1.0/f
        scale(i)=scale(i)*f
        noconv=.true.
        do 250 j=k,n
          ar(i,j)=ar(i,j)*g
          ai(i,j)=ai(i,j)*g
250     continue
        do 260 j=1,l
          ar(j,i)=ar(j,i)*f
          ai(j,i)=ai(j,i)*f
260     continue
270     continue
        if(noconv) goto 190
280     low=k
        igh=l
        return
        end

c
c    corth reduces a complex general matrix to upper hessenberg
c    form using unitary similarity transformations.
c
        subroutine corth(nm,n,low,igh,ar,ai,ortr,orti)
c
        integer i,j,m,n,ii,jj,la,mp,igh,nm,kp1,low
        real ar(nm,n),ai(nm,n),ortr(igh),orti(igh)
        real f,g,h,fi,fr,scale,sqrt,cabs,abs
        complex cmplx
c
        la=igh-1
        kp1=low+1
        if(la.lt.kp1) goto 200
        do 180 m=kp1,la
          h=0.0
          ortr(m)=0.0
          orti(m)=0.0
          scale=0.0
          do 90 i=m,igh
90              scale=scale+abs(ar(i,m-1))+abs(ai(i,m-1))
          if (scale.eq.0.0) goto 180
          mp=m+igh
          do 100 ii=m,igh
                i=mp-ii
                ortr(i)=ar(i,m-1)/scale
                orti(i)=ai(i,m-1)/scale
                h=h+ortr(i)*ortr(i)+orti(i)*orti(i)
100       continue
          g=sqrt(h)
          f=cabs(cmplx(ortr(m),orti(m)))
          if (f.eq.0.0) goto 103
          h=h+f*g
          g=g/f
          ortr(m)=(1.0+g)*ortr(m)
          orti(m)=(1.0+g)*orti(m)
          goto 105
103       ortr(m)=g
          ar(m,m-1)=scale
105       do 130 j=m,n
                fr=0.0
                fi=0.0
            do 110 ii=m,igh
                i=mp-ii
                fr=fr+ortr(i)*ar(i,j)+orti(i)*ai(i,j)
                fi=fi+ortr(i)*ai(i,j)-orti(i)*ar(i,j)
110         continue
                fr=fr/h
                fi=fi/h
            do 120 i=m,igh
                ar(i,j)=ar(i,j)-fr*ortr(i)+fi*orti(i)
                ai(i,j)=ai(i,j)-fr*orti(i)-fi*ortr(i)
120         continue
130       continue
          do 160 i=1,igh
                fr=0.0
                fi=0.0
            do 140 jj=m,igh
                j=mp-jj
                fr=fr+ortr(j)*ar(i,j)-orti(j)*ai(i,j)
                fi=fi+ortr(j)*ai(i,j)+orti(j)*ar(i,j)
140         continue
                fr=fr/h
                fi=fi/h
            do 150 j=m,igh
                ar(i,j)=ar(i,j)-fr*ortr(j)-fi*orti(j)
                ai(i,j)=ai(i,j)+fr*orti(j)-fi*ortr(j)
150         continue
160       continue
          ortr(m)=scale*ortr(m)
          orti(m)=scale*orti(m)
          ar(m,m-1)=-g*ar(m,m-1)
          ai(m,m-1)=-g*ai(m,m-1)
180     continue
200     return
        end

c
c   comqr2 use QR method to calculate eigenvalues and eigenvectors
c   from eispack.
c
      subroutine comqr2(nm,n,low,igh,ortr,orti,hr,hi,wr,wi,zr,zi,ierr)
        integer i,j,k,l,m,n,en,ii,jj,ll,nm,nn,igh,ip1,
     *    its,low,lp1,enm1,iend,ierr,min0
        real hr(nm,n),hi(nm,n),wr(n),wi(n),zr(nm,n),zi(nm,n),
     *    ortr(igh),orti(igh)
        real si,sr,ti,tr,xr,xi,yi,yr,zzi,zzr,norm,machep,abs
        complex z3,csqrt
c
c   machep, relative accuracy of single precision real number operation
c   on VAX 750
c
        machep=1.19e-07
        ierr=0
        do 100 i=1,n
          do 100 j=1,n
            zr(i,j)=0.0
            zi(i,j)=0.0
            if (i.eq.j) zr(i,j)=1.0
100     continue
        iend=igh-low-1
        if(iend) 180,150,105
105     do 140 ii=1,iend
          i=igh-ii
          if(ortr(i).eq.0.0.and.orti(i).eq.0.0) goto 140
          if(hr(i,i-1).eq.0.0.and.hi(i,i-1).eq.0.0) goto 140
          norm=hr(i,i-1)*ortr(i)+hi(i,i-1)*orti(i)
          ip1=i+1
          do 110 k=ip1,igh
                ortr(k)=hr(k,i-1)
                orti(k)=hi(k,i-1)
110       continue
          do 130 j=i,igh
                sr=0.0
                si=0.0
            do 115 k=i,igh
                sr=sr+ortr(k)*zr(k,j)+orti(k)*zi(k,j)
                si=si+ortr(k)*zi(k,j)-orti(k)*zr(k,j)
115         continue
                sr=sr/norm
                si=si/norm
            do 120 k=i,igh
                zr(k,j)=zr(k,j)+sr*ortr(k)-si*orti(k)
                zi(k,j)=zi(k,j)+sr*orti(k)+si*ortr(k)
120         continue
130       continue
140     continue
150     l=low+1
        do 170 i=l,igh
          ll=min0(i+1,igh)
          if(hi(i,i-1).eq.0.0) goto 170
          norm=cabs(cmplx(hr(i,i-1),hi(i,i-1)))
          yr=hr(i,i-1)/norm
          yi=hi(i,i-1)/norm
          hr(i,i-1)=norm
          hi(i,i-1)=0.0
          do 155 j=i,n
                si=yr*hi(i,j)-yi*hr(i,j)
                hr(i,j)=yr*hr(i,j)+yi*hi(i,j)
                hi(i,j)=si
155       continue
          do 160 j=1,ll
                si=yr*hi(j,i)+yi*hr(j,i)
                hr(j,i)=yr*hr(j,i)-yi*hi(j,i)
                hi(j,i)=si
160       continue
          do 165 j=low,igh
                si=yr*zi(j,i)+yi*zr(j,i)
                zr(j,i)=yr*zr(j,i)-yi*zi(j,i)
                zi(j,i)=si
165       continue
170     continue
180     do 200 i=1,n
          if (i.ge.low.and.i.le.igh) goto 200
          wr(i)=hr(i,i)
          wi(i)=hi(i,i)
200     continue
        en=igh
        tr=0.0
        ti=0.0
220     if (en.lt.low) goto 680
        its=0
        enm1=en-1
240     do 260 ll=low,en
          l=en+low-ll
          if (l.eq.low) goto 300
        if (abs(hr(l,l-1)).le.machep*(abs(hr(l-1,l-1))
     *  +abs(hi(l-1,l-1))+abs(hr(l,l))+abs(hi(l,l)))) goto 300
260     continue
300     if (l.eq.en) goto 660
        if (its.eq.30) goto 1000
        if (its.eq.10.or.its.eq.20) goto 320
        sr=hr(en,en)
        si=hi(en,en)
        xr=hr(enm1,en)*hr(en,enm1)
        xi=hi(enm1,en)*hr(en,enm1)
        if (xr.eq.0.0.and.xi.eq.0.0) goto 340
        yr=(hr(enm1,enm1)-sr)/2.0
        yi=(hi(enm1,enm1)-si)/2.0
        z3=csqrt(cmplx(yr**2-yi**2+xr,2.0*yr*yi+xi))
        zzr=real(z3)
        zzi=aimag(z3)
        if (yr*zzr+yi*zzi.ge.0.0) goto 310
        zzr=-zzr
        zzi=-zzi
310     z3=cmplx(xr,xi)/cmplx(yr+zzr,yi+zzi)
        sr=sr-real(z3)
        si=si-aimag(z3)
        goto 340
320     sr=abs(hr(en,enm1))+abs(hr(enm1,en-2))
        si=0.0
340     do 360 i=low,en
          hr(i,i)=hr(i,i)-sr
          hi(i,i)=hi(i,i)-si
360     continue
        tr=tr+sr
        ti=ti+si
        its=its+1
        lp1=l+1
        do 500 i=lp1,en
          sr=hr(i,i-1)
          hr(i,i-1)=0.0
          norm=sqrt(hr(i-1,i-1)*hr(i-1,i-1)+hi(i-1,i-1)
     *      *hi(i-1,i-1)+sr*sr)
          xr=hr(i-1,i-1)/norm
          wr(i-1)=xr
          xi=hi(i-1,i-1)/norm
          wi(i-1)=xi
          hr(i-1,i-1)=norm
          hi(i-1,i-1)=0.0
          hi(i,i-1)=sr/norm
c
          do 490 j=i,n
                yr=hr(i-1,j)
                yi=hi(i-1,j)
                zzr=hr(i,j)
                zzi=hi(i,j)
                hr(i-1,j)=xr*yr+xi*yi+hi(i,i-1)*zzr
                hi(i-1,j)=xr*yi-xi*yr+hi(i,i-1)*zzi
                hr(i,j)=xr*zzr-xi*zzi-hi(i,i-1)*yr
                hi(i,j)=xr*zzi+xi*zzr-hi(i,i-1)*yi
490       continue
500     continue
        si=hi(en,en)
        if(si.eq.0.0) goto 540
        norm=cabs(cmplx(hr(en,en),si))
        sr=hr(en,en)/norm
        si=si/norm
        hr(en,en)=norm
        hi(en,en)=0.0
        if(en.eq.n) goto 540
        ip1=en+1
        do 520 j=ip1,n
          yr=hr(en,j)
          yi=hi(en,j)
          hr(en,j)=sr*yr+si*yi
          hi(en,j)=sr*yi-si*yr
520     continue
540     do 600 j=lp1,en
          xr=wr(j-1)
          xi=wi(j-1)
          do 580 i=1,j
                yr=hr(i,j-1)
                yi=0.0
                zzr=hr(i,j)
                zzi=hi(i,j)
                if(i.eq.j) goto 560
                yi=hi(i,j-1)
                hi(i,j-1)=xr*yi+xi*yr+hi(j,j-1)*zzi
560             hr(i,j-1)=xr*yr-xi*yi+hi(j,j-1)*zzr
                hr(i,j)=xr*zzr+xi*zzi-hi(j,j-1)*yr
                hi(i,j)=xr*zzi-xi*zzr-hi(j,j-1)*yi
580       continue
          do 590 i=low,igh
                yr=zr(i,j-1)
                yi=zi(i,j-1)
                zzr=zr(i,j)
                zzi=zi(i,j)
                zr(i,j-1)=xr*yr-xi*yi+hi(j,j-1)*zzr
                zi(i,j-1)=xr*yi+xi*yr+hi(j,j-1)*zzi
                zr(i,j)=xr*zzr+xi*zzi-hi(j,j-1)*yr
                zi(i,j)=xr*zzi-xi*zzr-hi(j,j-1)*yi
590       continue
600     continue
        if(si.eq.0.0) goto 240
        do 630 i=1,en
          yr=hr(i,en)
          yi=hi(i,en)
          hr(i,en)=sr*yr-si*yi
          hi(i,en)=sr*yi+si*yr
630     continue
        do 640 i=low,igh
          yr=zr(i,en)
          yi=zi(i,en)
          zr(i,en)=sr*yr-si*yi
          zi(i,en)=sr*yi+si*yr
640     continue
        goto 240
660     hr(en,en)=hr(en,en)+tr
        wr(en)=hr(en,en)
        hi(en,en)=hi(en,en)+ti
        wi(en)=hi(en,en)
        en=enm1
        goto 220
680     norm=0.0
        do 720 i=1,n
          do 720 j=i,n
                norm=norm+abs(hr(i,j))+abs(hi(i,j))
720     continue
        if(n.eq.1.or.norm.eq.0.0) goto 1001
        do 800 nn=2,n
          en=n+2-nn
          xr=wr(en)
          xi=wi(en)
          enm1=en-1
          do 780 ii=1,enm1
                i=en-ii
                zzr=hr(i,en)
                zzi=hi(i,en)
                if(i.eq.enm1) goto 760
                ip1=i+1
            do 740 j=ip1,enm1
                zzr=zzr+hr(i,j)*hr(j,en)-hi(i,j)*hi(j,en)
                zzi=zzi+hr(i,j)*hi(j,en)+hi(i,j)*hr(j,en)
740         continue
760             yr=xr-wr(i)
                yi=xi-wi(i)
                if (yr.eq.0.0.and.yi.eq.0.0) yr=machep*norm
                z3=cmplx(zzr,zzi)/cmplx(yr,yi)
                hr(i,en)=real(z3)
                hi(i,en)=aimag(z3)
780       continue
800     continue
        enm1=n-1
        do 840 i=1,enm1
          if(i.ge.low.and.i.le.igh) goto 840
          ip1=i+1
          do 820 j=ip1,n
                zr(i,j)=hr(i,j)
                zi(i,j)=hi(i,j)
820       continue
840     continue
        do 880 jj=low,enm1
          j=n+low-jj
          m=min0(j-1,igh)
          do 880 i=low,igh
                zzr=zr(i,j)
                zzi=zi(i,j)
            do 860 k=low,m
                zzr=zzr+zr(i,k)*hr(k,j)-zi(i,k)*hi(k,j)
                zzi=zzi+zr(i,k)*hi(k,j)+zi(i,k)*hr(k,j)
860         continue
                zr(i,j)=zzr
                zi(i,j)=zzi
880     continue
        goto 1001
1000    ierr=en
1001    return
        end


        subroutine cbabk2(nm,n,low,igh,scale,m,zr,zi)
c
        real scale(n),zr(nm,m),zi(nm,m),s
        if (igh.eq.low) goto 120
        do 110 i=low,igh
          s=scale(i)
          do 100 j=1,m
                zr(i,j)=zr(i,j)*s
                zi(i,j)=zi(i,j)*s
100       continue
110     continue
120     do 140 ii=1,n
          i=ii
          if (i.ge.low.and.i.le.igh) goto 140
          if (i.lt.low) i=low-ii
          k=scale(i)
          if (k.eq.i) goto 140
          do 130 j=1,m
                s=zr(i,j)
                zr(i,j)=zr(k,j)
                zr(k,j)=s
                s=zi(i,j)
                zi(i,j)=zi(k,j)
                zi(k,j)=s
130       continue
140     continue
        return
        end
