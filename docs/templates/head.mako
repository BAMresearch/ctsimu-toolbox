<%!
    from pdoc.html_helpers import minify_css
%>
<%def name="homelink()" filter="minify_css">
    .homelink {
        display: block;
        font-size: 2em;
        font-weight: bold;
        color: #555;
        padding-bottom: .5em;
        border-bottom: 1px solid silver;
    }
    .homelink:hover {
        color: inherit;
    }
    .homelink img {
        max-width:20%;
        max-height: 5em;
        margin: auto;
        margin-bottom: .3em;
    }
</%def>

<style>${homelink()}</style>
<link rel="icon" href="https://bamresearch.github.io/ctsimu-toolbox/toolbox.png">
<style>
thead {
    border-top: 2px solid #000;
    border-bottom:1px solid #000;
}
tbody {
    border-bottom: 2px solid #000;
}
th, td {
    padding: 0.1em;
    padding-right: 1em;
    white-space: nowrap;
}
tr:hover {background-color: #eee;}
</style>