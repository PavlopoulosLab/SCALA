// This function disappears the requested button so the user does not
// click it multiple times before current function execution is finished
// @param btn: button element id
function disableButton(btn){
  var button = document.getElementById(btn);
  button.style.display = "none";
  return true;
}

// This function re-appears the requested button
// after the end of fucntion execution
// @param btn: button element id
function enableButton(btn){
  var button = document.getElementById(btn);
  button.style.display = "inline-block";
  return true;
}

// This function starts a requested loader
// @param m[0]: loader id
// @param m[1]: percentage to fill loader
function startLoader(m){
  var loader = document.getElementById(m[0]);
  loader.nextElementSibling.style.opacity = 0.5;
  loader.style.opacity = 1;
  loader.style.display = "inline-block";
  var bar = loader.ldBar;
  bar.set(parseInt(m[1]));
  return true;
}

// This function ends a requested loader
// @param m: loader id
function finishLoader(m){
  var loader = document.getElementById(m);
  var bar = loader.ldBar;
  bar.set(0);
  loader.nextElementSibling.style.opacity = 1;
  loader.style.display = "none";
  return true;
}

Shiny.addCustomMessageHandler("handler_disableButton", disableButton);
Shiny.addCustomMessageHandler("handler_enableButton", enableButton);
Shiny.addCustomMessageHandler("handler_startLoader", startLoader);
Shiny.addCustomMessageHandler("handler_finishLoader", finishLoader);
