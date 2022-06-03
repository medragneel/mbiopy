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


const rnaSeq = document.querySelector('#rnaSeq')
const protSeq = document.querySelector('#protSeq')
const mut = document.querySelector('#mutation')
const mut_count = document.querySelector('#mut_count')

async function saveRna() {
    fn=await window.showSaveFilePicker(); 
   let wr= await fn.createWritable();
    await wr.write(rnaSeq.innerText);
    await wr.close()
    
}


async function saveProt() {
    fn=await window.showSaveFilePicker(); 
   let wr= await fn.createWritable();
    await wr.write(protSeq.innerText);
    await wr.close()
    
}
async function saveMut() {
    fn=await window.showSaveFilePicker(); 
   let wr= await fn.createWritable();
    await wr.write(mut.innerText);
    await wr.write(`\n${mut_count.innerHTML}`);
    await wr.close()
    
}
