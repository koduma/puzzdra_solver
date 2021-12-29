function buttonClick(){
  msg.innerText = 'お名前は' + nameText.value + 'さんですね';
}

let nameText = document.getElementById('nameText');
nameText.value = '名前';
let msg = document.getElementById('msg');

let checkButton = document.getElementById('checkButton');
checkButton.addEventListener('click', butotnClick);
