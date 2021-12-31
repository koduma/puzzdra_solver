<!DOCTYPE html>
<html lang="ja">
<head>
<meta charset="UTF-8">
<title>パズドラ定石メーカーのルート計算サイト</title>
<script type="text/javascript" src="https://ajax.googleapis.com/ajax/libs/jquery/3.6.0/jquery.min.js"></script>
<script type="text/javascript">
$('form').submit(function(){

    const myObj = {code: '#include <iostream>\nint main() { int x = 0; std::cout << \"hoge\" << std::endl; }',options:'warning,gnu++1y',
                   compiler:'gcc-head',compiler-option-raw:'-Dx=hogefuga\n-O3'};

    $.ajax({
        type: 'POST',
        url: 'https://wandbox.org/api/compile.json',
        dataType: 'json',
        data: JSON.stringify(myObj),
        contentType : 'application/json',
        success: function(data) {
        successmessage = 'Data was succesfully captured';
        $("label#successmessage").text(successmessage);
        },
        error: function (data) {
        successmessage = 'Error';
        $("label#successmessage").text(successmessage);
        },
    });
});
</script>

</head>
<body>

<p>パズドラ定石メーカーのルート計算サイトです。kkk</p>

<form>
<label>layout=<input type="text" id="layout"></label>
<button value="submit" type="submit">計算</button>
</form>
  
<p id="msg"></p>
  
<a href="#" onclick="redirect();">計算結果</a>

<script>
  
function redirect(){
    location.href = 'http://serizawa.web5.jp/puzzdra_theory_maker/index.html?layout=' + layout.value + '&route=null&ctwMode=false';
}

let layout = document.getElementById('layout');
layout.value = '';
let msg = document.getElementById('msg');
let checkButton = document.getElementById('checkButton');
checkButton.addEventListener('click', buttonClick);
</script>

</body>
</html>
