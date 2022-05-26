window.addEventListener('scroll',function(){
    const nav = document.querySelector('nav');
    nav.classList.toggle('on_scroll',window.scrollY > 0);
    document.getElementById('logo').src="../static/images/logo-dark-bio-py-tech.png";
                })


bg= document.querySelector('.ham')
ul= document.querySelector('.nav-links')
bg.addEventListener('click',()=>{
    ul.classList.toggle('active')
                })

