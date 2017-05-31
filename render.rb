#!/usr/bin/env ruby

# sudo apt-get install ruby-dev
# sudo gem install github-markdown

# render.rb
require 'github/markdown'

body = GitHub::Markdown.render_gfm File.read(ARGV[0])

file = File.open("github-markdown.css", "rb")
css = file.read

puts '<html>
<head>
<title>Replication Package for Very Simple Markov-Perfect Industry Dynamics: Theory </title>
<style>' +
css + '
	.markdown-body {
		box-sizing: border-box;
		min-width: 200px;
		max-width: 980px;
		margin: 0 auto;
		padding: 45px;
	}
</style>
</head>
<body>
<article class="markdown-body">' + body + '</article>
</body>
</html>'
